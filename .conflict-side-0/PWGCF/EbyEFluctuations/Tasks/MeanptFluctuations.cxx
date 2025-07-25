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

/// \file MeanptFluctuations.cxx
/// \brief Task for analyzing <pT> fluctuation upto fourth order of inclusive hadrons
/// \author Swati Saha

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <string>
#include <vector>

namespace o2::aod
{
namespace ptQn
{
DECLARE_SOA_COLUMN(Q1, q1, float);                 //! sum of pT of tracks in an event
DECLARE_SOA_COLUMN(Q2, q2, float);                 //! sum of (pT)^2 of tracks in an event
DECLARE_SOA_COLUMN(Q3, q3, float);                 //! sum of (pT)^3 of tracks in an event
DECLARE_SOA_COLUMN(Q4, q4, float);                 //! sum of (pT)^4 of tracks in an event
DECLARE_SOA_COLUMN(N_ch, n_ch, float);             //! no of charged particles/multiplicity in an event
DECLARE_SOA_COLUMN(Centrality, centrality, float); //! Centrality of event
} // namespace ptQn
DECLARE_SOA_TABLE(MultPtQn, "AOD", "PTQN", ptQn::Q1, ptQn::Q2, ptQn::Q3, ptQn::Q4, ptQn::N_ch, ptQn::Centrality); //! table to store e-by-e sum of pT, (pT)^2, (pT)^3, (pT)^4 of tracks, multiplicity and centrality
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct MeanptFluctuations_QA_QnTable {

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 3.0f, "Higher pT cut"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
  Configurable<float> cfgCutItsChi2NCl{"cfgCutItsChi2NCl", 36.0f, "Maximum ITSchi2NCl"};
  Configurable<float> cfgCutTrackDcaZ{"cfgCutTrackDcaZ", 2.0f, "Maximum DcaZ"};
  Configurable<int> cfgITScluster{"cfgITScluster", 1, "Minimum Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 80, "Minimum Number of TPC cluster"};
  Configurable<int> cfgTPCnCrossedRows{"cfgTPCnCrossedRows", 70, "Minimum Number of TPC crossed-rows"};
  ConfigurableAxis nchAxis{"nchAxis", {5000, 0.5, 5000.5}, ""};
  Configurable<bool> cfgEvSelkNoSameBunchPileup{"cfgEvSelkNoSameBunchPileup", true, "Pileup removal"};
  Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Remove time interval with dead ITS zone"};

  O2_DEFINE_CONFIGURABLE(cfgUse22sEventCut, bool, true, "Use 22s event cut on mult correlations")

  // Filter command***********
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < 0.8f) && (aod::track::pt > cfgCutPtLower) && (aod::track::pt < 5.0f) && (requireGlobalTrackInFilter()) && (aod::track::tpcChi2NCl < cfgCutTpcChi2NCl) && (aod::track::itsChi2NCl < cfgCutItsChi2NCl) && (nabs(aod::track::dcaZ) < cfgCutTrackDcaZ);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // filtering collisions and tracks***********
  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>>;
  // using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    std::vector<double> centBining = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    AxisSpec centAxis = {centBining, "centrality (%)"};

    // Add histograms to histogram manager (as in the output object of in AliPhysics)
    histos.add("hZvtx_after_sel", ";Z (cm)", kTH1F, {vtxZAxis});
    histos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.2, 4.}});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hPhi", ";#phi", kTH1F, {{100, 0., o2::constants::math::TwoPI}});
    histos.add("hEta", ";#eta", kTH1F, {{100, -2.01, 2.01}});
    histos.add("hCentrality", ";centrality (%)", kTH1F, {{90, 0, 90}});
    histos.add("hDcaXY", ";#it{dca}_{XY}", kTH1F, {{1000, -5, 5}});
    histos.add("hDcaZ", ";#it{dca}_{Z}", kTH1F, {{1000, -5, 5}});
    histos.add("hMeanPt", "", kTProfile, {centAxis});
    histos.add("Hist2D_globalTracks_PVTracks", "", {HistType::kTH2D, {nchAxis, nchAxis}});
    histos.add("Hist2D_cent_nch", "", {HistType::kTH2D, {nchAxis, centAxis}});

    // Event selection - Alex
    if (cfgUse22sEventCut) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
      fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
      fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
      fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
    }

  } //! end init function

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk, const float& centrality)
  {
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      return 0;
    }
    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      if (zRes > 0.25 && collision.numContrib() < 20)
        vtxz = -999;
    }
    auto multNTracksPV = collision.multNTracksPV();

    if ((vtxz > cfgCutVertex) || (vtxz < -1.0 * cfgCutVertex))
      return 0;
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;
    if (multTrk < fMultCutLow->Eval(centrality))
      return 0;
    if (multTrk > fMultCutHigh->Eval(centrality))
      return 0;
    if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
      return 0;

    return 1;
  }

  Produces<aod::MultPtQn> mult_ptQn;

  // void process(aod::Collision const& coll, aod::Tracks const& inputTracks)
  void process(aodCollisions::iterator const& coll, aod::BCsWithTimestamps const&, aodTracks const& inputTracks)
  {
    if (!coll.sel8()) {
      return;
    }
    if (cfgUseGoodITSLayerAllCut && !(coll.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))) {
      return;
    }
    if (cfgEvSelkNoSameBunchPileup && !(coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
      return;
    }

    const auto CentralityFT0C = coll.centFT0C();
    if (cfgUse22sEventCut && !eventSelected(coll, inputTracks.size(), CentralityFT0C))
      return;

    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), coll.centFT0C());
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), inputTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), inputTracks.size(), CentralityFT0C);

    // variables
    double cent = coll.centFT0C();
    double pT_sum = 0.0;
    double N = 0.0;

    float q1 = 0.0;
    float q2 = 0.0;
    float q3 = 0.0;
    float q4 = 0.0;
    float n_ch = 0.0;

    for (const auto& track : inputTracks) { // Loop over tracks

      if (!track.has_collision()) {
        continue;
      }

      if (!track.isPVContributor()) {
        continue;
      }

      if (!(track.itsNCls() > cfgITScluster) || !(track.tpcNClsFound() >= cfgTPCcluster) || !(track.tpcNClsCrossedRows() >= cfgTPCnCrossedRows)) {
        continue;
      }

      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hPhi"), track.phi());
      histos.fill(HIST("hDcaXY"), track.dcaXY());
      histos.fill(HIST("hDcaZ"), track.dcaZ());

      pT_sum += track.pt();
      N += 1.0;

      float pT = track.pt();
      // calculating Q1, Q2, Q3, Q4. N_ch
      if (track.pt() > cfgCutPtLower && track.pt() < cfgCutPtUpper && track.sign() != 0) {
        q1 = q1 + std::pow(pT, 1.0);
        q2 = q2 + std::pow(pT, 2.0);
        q3 = q3 + std::pow(pT, 3.0);
        q4 = q4 + std::pow(pT, 4.0);
        n_ch = n_ch + 1;
      }
    }
    mult_ptQn(q1, q2, q3, q4, n_ch, cent);
    // MeanPt
    if (N > 0.0f)
      histos.fill(HIST("hMeanPt"), cent, pT_sum / N);
  }
};

struct MeanptFluctuations_analysis {

  Configurable<int> cfgNSubsample{"cfgNSubsample", 10, "Number of subsamples"};
  ConfigurableAxis centAxis{"centAxis", {90, 0, 90}, ""};
  ConfigurableAxis multAxis{"multAxis", {5000, 0.5, 5000.5}, ""};
  ConfigurableAxis meanpTAxis{"meanpTAxis", {500, 0, 5.0}, ""};

  expressions::Filter Nch_filter = aod::ptQn::n_ch > 3.0f;
  using FilteredMultPtQn = soa::Filtered<aod::MultPtQn>;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // Define output
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<std::vector<std::shared_ptr<TProfile2D>>> Subsample;
  TRandom3* fRndm = new TRandom3(0);

  void init(o2::framework::InitContext&)
  {
    // AxisSpec centAxis = {90, 0, 90, "centrality (%)"};
    // AxisSpec multAxis = {5000, 0.5, 5000.5, "#it{N}_{ch,acc}"};

    registry.add("Prof_mean_t1", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof_var_t1", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof_skew_t1", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof_kurt_t1", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Hist2D_Nch_centrality", "", {HistType::kTH2D, {centAxis, multAxis}});
    registry.add("Hist2D_meanpt_centrality", "", {HistType::kTH2D, {centAxis, meanpTAxis}});

    // initial array
    Subsample.resize(cfgNSubsample);
    for (int i = 0; i < cfgNSubsample; i++) {
      Subsample[i].resize(4);
    }
    for (int i = 0; i < cfgNSubsample; i++) {
      Subsample[i][0] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof_mean_t1", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      Subsample[i][1] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof_var_t1", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      Subsample[i][2] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof_skew_t1", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      Subsample[i][3] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof_kurt_t1", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
    }
  }

  float mean_term1;
  float variance_term1;
  float skewness_term1;
  float kurtosis_term1;

  // void process(aod::MultPtQn::iterator const& event_ptqn)
  void process(FilteredMultPtQn::iterator const& event_ptqn)
  {
    // LOGF(info, "Centrality= %f Nch= %f Q1= %f Q2= %f", event_ptqn.centrality(), event_ptqn.n_ch(), event_ptqn.q1(), event_ptqn.q2());

    // calculating observables
    mean_term1 = event_ptqn.q1() / event_ptqn.n_ch();
    variance_term1 = (std::pow(event_ptqn.q1(), 2.0f) - event_ptqn.q2()) / (event_ptqn.n_ch() * (event_ptqn.n_ch() - 1.0f));
    skewness_term1 = (std::pow(event_ptqn.q1(), 3.0f) - 3.0f * event_ptqn.q2() * event_ptqn.q1() + 2.0f * event_ptqn.q3()) / (event_ptqn.n_ch() * (event_ptqn.n_ch() - 1.0f) * (event_ptqn.n_ch() - 2.0f));
    kurtosis_term1 = (std::pow(event_ptqn.q1(), 4.0f) - (6.0f * event_ptqn.q4()) + (8.0f * event_ptqn.q1() * event_ptqn.q3()) - (6.0f * std::pow(event_ptqn.q1(), 2.0f) * event_ptqn.q2()) + (3.0f * std::pow(event_ptqn.q2(), 2.0f))) / (event_ptqn.n_ch() * (event_ptqn.n_ch() - 1.0f) * (event_ptqn.n_ch() - 2.0f) * (event_ptqn.n_ch() - 3.0f));

    // filling profiles and histograms for central values
    registry.get<TProfile2D>(HIST("Prof_mean_t1"))->Fill(event_ptqn.centrality(), event_ptqn.n_ch(), mean_term1);
    registry.get<TProfile2D>(HIST("Prof_var_t1"))->Fill(event_ptqn.centrality(), event_ptqn.n_ch(), variance_term1);
    registry.get<TProfile2D>(HIST("Prof_skew_t1"))->Fill(event_ptqn.centrality(), event_ptqn.n_ch(), skewness_term1);
    registry.get<TProfile2D>(HIST("Prof_kurt_t1"))->Fill(event_ptqn.centrality(), event_ptqn.n_ch(), kurtosis_term1);
    registry.fill(HIST("Hist2D_Nch_centrality"), event_ptqn.centrality(), event_ptqn.n_ch());
    registry.fill(HIST("Hist2D_meanpt_centrality"), event_ptqn.centrality(), mean_term1);

    // selecting subsample and filling profiles
    float l_Random = fRndm->Rndm();
    int SampleIndex = static_cast<int>(cfgNSubsample * l_Random);
    Subsample[SampleIndex][0]->Fill(event_ptqn.centrality(), event_ptqn.n_ch(), mean_term1);
    Subsample[SampleIndex][1]->Fill(event_ptqn.centrality(), event_ptqn.n_ch(), variance_term1);
    Subsample[SampleIndex][2]->Fill(event_ptqn.centrality(), event_ptqn.n_ch(), skewness_term1);
    Subsample[SampleIndex][3]->Fill(event_ptqn.centrality(), event_ptqn.n_ch(), kurtosis_term1);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  return WorkflowSpec{
    adaptAnalysisTask<MeanptFluctuations_QA_QnTable>(cfgc),
    adaptAnalysisTask<MeanptFluctuations_analysis>(cfgc),
  };
}
