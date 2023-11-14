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

// Task to produce a table joinable to the jet tables for hf jet tagging
//
// Author: Hanseo Park

#include "TF1.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/JetTaggingUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTagTable>
struct JetTaggerHFTC {

  static constexpr std::string_view charJetFlavor[] = {"inc_jet", "lfjet", "cjet", "bjet"};

  Configurable<bool> doFitResoFunc{"doFitResoFunc", false, "Do Fit resolution function"};
  // Binning
  ConfigurableAxis binJetPt{"binJetPt", {200, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binNtracks{"binNtracks", {100, -0.5, 99.5}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binImpactParameterXY{"binImpactParameterXY", {1000, -0.4f, 0.4f}, ""};
  ConfigurableAxis binImpactParameterSignificanceXY{"binImpactParameterSignificanceXY", {1000, -40.f, 40.f}, ""};
  ConfigurableAxis binImpactParameterXYZ{"binImpactParameterXYZ", {1000, -0.4f, 0.4f}, ""};
  ConfigurableAxis binImpactParameterSignificanceXYZ{"binImpactParameterSignificanceXYZ", {1000, -40.f, 40.f}, ""};
  ConfigurableAxis binTrackProbability{"binTrackProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetProbability{"binJetProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetProbabilityLog{"binJetProbabilityLog", {100, 0.f, 10.f}, ""};
  ConfigurableAxis binJetEntries{"binEntries", {3, 0.f, 3.f}, ""};
  Configurable<std::vector<double>> tempbinsJetPt{"tempbinsJetPt", std::vector<double>{JetTaggingBinCut::vecBinsJetPt}, "pT bin limits"};

  // Axis
  AxisSpec jetPtRangeAxis{(std::vector<double>)tempbinsJetPt, "#it{p}_{T}^{Jet} (GeV/#it{c})"};
  AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T, jet}"};
  AxisSpec etaAxis = {binEta, "#eta"};
  AxisSpec phiAxis = {binPhi, "#phi"};
  AxisSpec ntracksAxis = {binNtracks, "#it{N}_{tracks}"};
  AxisSpec trackPtAxis = {binTrackPt, "#it{p}_{T}^{track}"};
  AxisSpec ImpactParameterXYAxis = {binImpactParameterXY, "IP_{XY} [cm]"};
  AxisSpec ImpactParameterSignificanceXYAxis = {binImpactParameterSignificanceXY, "IPs_{XY} [cm]"};
  AxisSpec ImpactParameterXYZAxis = {binImpactParameterXYZ, "IP_{XYZ} [cm]"};
  AxisSpec ImpactParameterSignificanceXYZAxis = {binImpactParameterSignificanceXYZ, "IPs_{XYZ} [cm]"};
  AxisSpec TrackProbabilityAxis = {binTrackProbability, "track probability"};
  AxisSpec JetProbabilityAxis = {binJetProbability, "JP"};
  AxisSpec JetProbabilityLogAxis = {binJetProbabilityLog, "-Log(JP)"};
  AxisSpec JetEntries = {binJetEntries, "lf=1, c=2, b=3"};

  TF1* fResoFunccjet = new TF1("fResoFunccjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  TF1* fResoFuncbjet = new TF1("fResoFuncbjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  TF1* fResoFunclfjet = new TF1("fResoFunclfjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  TF1* fResoFuncincjet = new TF1("fResoFuncincjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    for (int i = 0; i < 4; i++) {
      // common
      registry.add(Form("h_%s_pt_jet_eta", charJetFlavor[i].data()), Form("%s #it{p}_{T}", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, etaAxis}}, true);
      registry.add(Form("h_%s_pt", charJetFlavor[i].data()), Form("%s #it{p}_{T}", charJetFlavor[i].data()), {HistType::kTH1F, {jetPtAxis}}, true);
      registry.add(Form("h_%s_eta", charJetFlavor[i].data()), Form("%s #eta", charJetFlavor[i].data()), {HistType::kTH1F, {etaAxis}}, true);
      registry.add(Form("h_%s_phi", charJetFlavor[i].data()), Form("%s #phi", charJetFlavor[i].data()), {HistType::kTH1F, {phiAxis}}, true);
      registry.add(Form("h_%s_ntracks", charJetFlavor[i].data()), Form("%s N tracks", charJetFlavor[i].data()), {HistType::kTH1F, {ntracksAxis}}, true);
      registry.add(Form("h_%s_track_pt", charJetFlavor[i].data()), Form("#it{p}_{T} of track in %s", charJetFlavor[i].data()), {HistType::kTH1F, {trackPtAxis}}, true);
      registry.add(Form("h_%s_track_eta", charJetFlavor[i].data()), Form("#eta of track in %s", charJetFlavor[i].data()), {HistType::kTH1F, {etaAxis}}, true);
      registry.add(Form("h_%s_track_phi", charJetFlavor[i].data()), Form("#phi of track in %s", charJetFlavor[i].data()), {HistType::kTH1F, {phiAxis}}, true);
      registry.add(Form("h2_%s_track_pt", charJetFlavor[i].data()), Form("#it{p}_{T} of track in %s", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtAxis, trackPtAxis}}, true);
      registry.add(Form("h2_%s_track_eta", charJetFlavor[i].data()), Form("#it{p}_{T} of track in %s", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtAxis, etaAxis}}, true);
      registry.add(Form("h2_%s_track_phi", charJetFlavor[i].data()), Form("#it{p}_{T} of track in %s", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtAxis, phiAxis}}, true);
      registry.add(Form("h_%s_impact_parameter_xy", charJetFlavor[i].data()), Form("%s impact parameter dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xy", charJetFlavor[i].data()), Form("%s sign impact parameter dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYAxis}});
      registry.add(Form("h_%s_impact_parameter_xy_significance", charJetFlavor[i].data()), Form("%s impact parameter significance dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xy_significance", charJetFlavor[i].data()), Form("%s sign impact parameter significance dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
      registry.add(Form("h_%s_pt_sign_impact_parameter_xy_significance", charJetFlavor[i].data()), Form("%s sign impact parameter significance dca_{xy}", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, ImpactParameterSignificanceXYAxis}});
      registry.add(Form("h_%s_impact_parameter_xyz", charJetFlavor[i].data()), Form("%s impact parameter dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYZAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xyz", charJetFlavor[i].data()), Form("%s sign impact parameter dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYZAxis}});
      registry.add(Form("h_%s_impact_parameter_xyz_significance", charJetFlavor[i].data()), Form("%s impact parameter significance dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xyz_significance", charJetFlavor[i].data()), Form("%s sign impact parameter significance dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
      registry.add(Form("h_%s_pt_sign_impact_parameter_xyz_significance", charJetFlavor[i].data()), Form("%s sign impact parameter significance dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, ImpactParameterSignificanceXYZAxis}});

      for (int j = 1; j < 4; j++) { // Track counting
        registry.add(Form("h_%s_impact_parameter_xy_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xy_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYAxis}});
        registry.add(Form("h_%s_impact_parameter_xy_significance_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter significance dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xy_significance_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
        registry.add(Form("h_%s_pt_sign_impact_parameter_xy_significance_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH2F, {jetPtRangeAxis, ImpactParameterSignificanceXYAxis}});
        registry.add(Form("h_%s_impact_parameter_xyz_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYZAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xyz_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYZAxis}});
        registry.add(Form("h_%s_impact_parameter_xyz_significance_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter significance dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xyz_significance_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
        registry.add(Form("h_%s_pt_sign_impact_parameter_xyz_significance_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH2F, {jetPtRangeAxis, ImpactParameterSignificanceXYZAxis}});
      }
      if (doprocessJPMCD) {
        registry.add(Form("h_%s_pos_track_probability", charJetFlavor[i].data()), Form("%s positive sign track probability", charJetFlavor[i].data()), {HistType::kTH1F, {TrackProbabilityAxis}});
        registry.add(Form("h_%s_neg_track_probability", charJetFlavor[i].data()), Form("%s negative sign track probability", charJetFlavor[i].data()), {HistType::kTH1F, {TrackProbabilityAxis}});
        registry.add(Form("h_%s_pt_pos_track_probability", charJetFlavor[i].data()), Form("%s positive sign track probability", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, TrackProbabilityAxis}});
        registry.add(Form("h_%s_pt_neg_track_probability", charJetFlavor[i].data()), Form("%s negative sign track probability", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, TrackProbabilityAxis}});
        registry.add(Form("h_%s_JP", charJetFlavor[i].data()), Form("%s jet probability", charJetFlavor[i].data()), {HistType::kTH1F, {JetProbabilityAxis}});
        registry.add(Form("h_%s_JP_Log", charJetFlavor[i].data()), Form("Log %s jet probability", charJetFlavor[i].data()), {HistType::kTH1F, {JetProbabilityLogAxis}});
        registry.add(Form("h_%s_pt_JP", charJetFlavor[i].data()), Form("%s jet probability", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, JetProbabilityAxis}});
        registry.add(Form("h_%s_pt_JP_Log", charJetFlavor[i].data()), Form("Log %s jet probability", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, JetProbabilityLogAxis}});
      }
    }
  }

  using JetTagTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::StoredTracksIU, aod::TracksCovIU, aod::McTrackLabels>;

  void SetFitResoFunc()
  {
    if (doFitResoFunc) {
      // TODO
      LOGF(info, "Bring resolution function from CCDB");
    } else {
      fResoFuncincjet->SetParameters(1.28558e+00, 1.51891e-01, 3.27502e+00, 4.41780e-01, 6.08111e+00, 7.14751e-01, 4.31734e+03, 8.06120e-02, 9.71495e-01);  // refer sign to IPs distribution
      fResoFunclfjet->SetParameters(3.36843e+00, 3.96023e-01, 4.25250e-01, 7.29701e-02, 6.50267e+00, 8.82015e-01, 3.96449e+03, 3.45753e-02, 9.52941e-01);   // refer sign to IPs distribution
      fResoFunccjet->SetParameters(8.19589e-01, 6.82005e-02, 2.27212e+00, 1.69561e+00, 5.68792e+00, 1.63553e+00, -1.10237e+03, 1.43584e+00, 7.90208e-01);   // refer to sign IPs distribution
      fResoFuncbjet->SetParameters(-1.19101e+01, 2.48039e+01, 3.14354e-01, 3.12194e-02, 3.19051e+00, 1.22309e+00, -1.63482e+03, 1.24152e+00, -9.46069e-02); // refer to sign IPs distribution
    }
  }

  template <typename T>
  bool trackSelectionJP(T const& jet)
  {
    if (jet.tracks().size() < 2)
      return 0;

    return 1;
  }

  template <typename T, typename U, typename V, typename D = double>
  double CalculateJP(T const& collision, U const& jet, V const& mcParticles, D& jetflavor)
  {
    double JP = -1.;
    if (!trackSelectionJP(jet))
      return -1.;
    SetFitResoFunc();
    std::vector<double> jetTracksPt;

    double trackjetProb = 1.;
    if (jetflavor == JetTaggingSpecies::none)
      return -1;
    for (auto& track : jet.template tracks_as<JetTagTracks>()) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for track, skip...");
        continue;
      }
      auto particle = track.mcParticle();
      if (!particle.has_mothers()) {
        LOGF(warning, "No mother particle for particle, skip...");
        continue;
      }

      int sgn;
      double varImpXY, varSgnImpXY, varImpXYSig, varSgnImpXYSig, varImpXYZ, varSgnImpXYZ, varImpXYZSig, varSgnImpXYZSig;
      JetTaggingUtilities::SetSgnImpactParameterSignificance(collision, jet, track, sgn, varImpXY, varSgnImpXY, varImpXYSig, varSgnImpXYSig, varImpXYZ, varSgnImpXYZ, varImpXYZSig, varSgnImpXYZSig);

      if (jetflavor == JetTaggingSpecies::none) {
        LOGF(debug, "NOT DEFINE JET FLAVOR");
        continue;
      }

      double ProbTrack = 0.;
      auto minSgnImpXYSig = -20;
      ProbTrack = fResoFuncincjet->Integral(minSgnImpXYSig, -1 * TMath::Abs(varSgnImpXYSig)) / fResoFuncincjet->Integral(minSgnImpXYSig, 0);
      if (TMath::Abs(varSgnImpXYSig) > TMath::Abs(minSgnImpXYSig))
        varSgnImpXYSig = 9.99; // Limit to function definition range
      if (sgn > 0) {
        registry.fill(HIST("h_inc_jet_pos_track_probability"), ProbTrack);
        registry.fill(HIST("h_inc_jet_pt_pos_track_probability"), jet.pt(), ProbTrack);
      } else {
        registry.fill(HIST("h_inc_jet_neg_track_probability"), ProbTrack);
        registry.fill(HIST("h_inc_jet_pt_neg_track_probability"), jet.pt(), ProbTrack);
      }
      if (jetflavor == JetTaggingSpecies::lightflavour) {
        ProbTrack = fResoFunclfjet->Integral(minSgnImpXYSig, -1 * TMath::Abs(varSgnImpXYSig)) / fResoFunclfjet->Integral(minSgnImpXYSig, 0);
        if (sgn > 0) {
          registry.fill(HIST("h_lfjet_pos_track_probability"), ProbTrack);
          registry.fill(HIST("h_lfjet_pt_pos_track_probability"), jet.pt(), ProbTrack);
        } else {
          registry.fill(HIST("h_lfjet_neg_track_probability"), ProbTrack);
          registry.fill(HIST("h_lfjet_pt_neg_track_probability"), jet.pt(), ProbTrack);
        }
      }
      if (jetflavor == JetTaggingSpecies::charm) {
        ProbTrack = fResoFunccjet->Integral(minSgnImpXYSig, -1 * TMath::Abs(varSgnImpXYSig)) / fResoFunccjet->Integral(minSgnImpXYSig, 0);
        if (sgn > 0) {
          registry.fill(HIST("h_cjet_pos_track_probability"), ProbTrack);
          registry.fill(HIST("h_cjet_pt_pos_track_probability"), jet.pt(), ProbTrack);
        } else {
          registry.fill(HIST("h_cjet_neg_track_probability"), ProbTrack);
          registry.fill(HIST("h_cjet_pt_neg_track_probability"), jet.pt(), ProbTrack);
        }
      }
      if (jetflavor == JetTaggingSpecies::beauty) {
        ProbTrack = fResoFuncbjet->Integral(minSgnImpXYSig, -1 * TMath::Abs(varSgnImpXYSig)) / fResoFuncbjet->Integral(minSgnImpXYSig, 0);
        if (sgn > 0) {
          registry.fill(HIST("h_bjet_pos_track_probability"), ProbTrack);
          registry.fill(HIST("h_bjet_pt_pos_track_probability"), jet.pt(), ProbTrack);
        } else {
          registry.fill(HIST("h_bjet_neg_track_probability"), ProbTrack);
          registry.fill(HIST("h_bjet_pt_neg_track_probability"), jet.pt(), ProbTrack);
        }
      }
      if (sgn > 0) { // only take postive sign track
        trackjetProb = trackjetProb * TMath::Abs(ProbTrack);
        jetTracksPt.push_back(track.pt());
      }
    }
    double sumjetProb = 0.;

    if (jetTracksPt.size() < 2)
      return -1;

    for (int i = 0; i < jetTracksPt.size(); i++) { // JP
      sumjetProb = sumjetProb + (TMath::Power(-1 * TMath::Log(trackjetProb), i) / TMath::Factorial(i));
    }

    JP = trackjetProb * sumjetProb;

    return JP;
  }

  template <typename T, typename U>
  void fillHistogramIPsMCD(T const& collision, U const& mcdjets)
  {
    for (auto& mcdjet : mcdjets) {
      const int jetflavor = mcdjet.origin();
      if (jetflavor == JetTaggingSpecies::none) {
        LOGF(debug, "NOT DEFINE JET FLAVOR");
      }
      // common var
      std::vector<double> incjetTracksImpXY, lfjetTracksImpXY, cjetTracksImpXY, bjetTracksImpXY;
      std::vector<double> incjetTracksSignImpXY, lfjetTracksSignImpXY, cjetTracksSignImpXY, bjetTracksSignImpXY;
      std::vector<double> incjetTracksImpXYSig, lfjetTracksImpXYSig, cjetTracksImpXYSig, bjetTracksImpXYSig;
      std::vector<double> incjetTracksSignImpXYSig, lfjetTracksSignImpXYSig, cjetTracksSignImpXYSig, bjetTracksSignImpXYSig;
      std::vector<double> incjetTracksImpXYZ, lfjetTracksImpXYZ, cjetTracksImpXYZ, bjetTracksImpXYZ;
      std::vector<double> incjetTracksSignImpXYZ, lfjetTracksSignImpXYZ, cjetTracksSignImpXYZ, bjetTracksSignImpXYZ;
      std::vector<double> incjetTracksImpXYZSig, lfjetTracksImpXYZSig, cjetTracksImpXYZSig, bjetTracksImpXYZSig;
      std::vector<double> incjetTracksSignImpXYZSig, lfjetTracksSignImpXYZSig, cjetTracksSignImpXYZSig, bjetTracksSignImpXYZSig;

      registry.fill(HIST("h_inc_jet_pt_jet_eta"), mcdjet.pt(), mcdjet.eta());
      registry.fill(HIST("h_inc_jet_pt"), mcdjet.pt());
      registry.fill(HIST("h_inc_jet_eta"), mcdjet.eta());
      registry.fill(HIST("h_inc_jet_phi"), mcdjet.phi());
      registry.fill(HIST("h_inc_jet_ntracks"), mcdjet.tracks().size());

      if (jetflavor == JetTaggingSpecies::lightflavour) {
        registry.fill(HIST("h_lfjet_pt"), mcdjet.pt());
        registry.fill(HIST("h_lfjet_eta"), mcdjet.eta());
        registry.fill(HIST("h_lfjet_phi"), mcdjet.phi());
        registry.fill(HIST("h_lfjet_ntracks"), mcdjet.tracks().size());
      }
      if (jetflavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h_cjet_pt"), mcdjet.pt());
        registry.fill(HIST("h_cjet_eta"), mcdjet.eta());
        registry.fill(HIST("h_cjet_phi"), mcdjet.phi());
        registry.fill(HIST("h_cjet_ntracks"), mcdjet.tracks().size());
      }
      if (jetflavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_bjet_pt"), mcdjet.pt());
        registry.fill(HIST("h_bjet_eta"), mcdjet.eta());
        registry.fill(HIST("h_bjet_phi"), mcdjet.phi());
        registry.fill(HIST("h_bjet_ntracks"), mcdjet.tracks().size());
      }

      for (auto& track : mcdjet.template tracks_as<JetTagTracks>()) {
        int sign;
        double varImpXY, varSgnImpXY, varImpXYSig, varSgnImpXYSig, varImpXYZ, varSgnImpXYZ, varImpXYZSig, varSgnImpXYZSig;
        JetTaggingUtilities::SetSgnImpactParameterSignificance(collision, mcdjet, track, sign, varImpXY, varSgnImpXY, varImpXYSig, varSgnImpXYSig, varImpXYZ, varSgnImpXYZ, varImpXYZSig, varSgnImpXYZSig);
        registry.fill(HIST("h_inc_jet_track_pt"), track.pt());
        registry.fill(HIST("h_inc_jet_track_phi"), track.phi());
        registry.fill(HIST("h_inc_jet_track_eta"), track.eta());
        registry.fill(HIST("h_inc_jet_impact_parameter_xy"), varImpXY);
        registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy"), varSgnImpXY);
        registry.fill(HIST("h_inc_jet_impact_parameter_xy_significance"), varImpXYSig);
        registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_significance"), varSgnImpXYSig);
        registry.fill(HIST("h_inc_jet_pt_sign_impact_parameter_xy_significance"), mcdjet.pt(), varSgnImpXYSig);
        registry.fill(HIST("h_inc_jet_impact_parameter_xyz"), varImpXYZ);
        registry.fill(HIST("h_inc_jet_sign_impact_parameter_xyz"), varSgnImpXYZ);
        registry.fill(HIST("h_inc_jet_impact_parameter_xyz_significance"), varImpXYZSig);
        registry.fill(HIST("h_inc_jet_sign_impact_parameter_xyz_significance"), varSgnImpXYZSig);
        incjetTracksImpXY.push_back(varImpXY);
        incjetTracksImpXYSig.push_back(varImpXYSig);
        incjetTracksSignImpXY.push_back(varSgnImpXY);
        incjetTracksSignImpXYSig.push_back(varSgnImpXYSig);
        incjetTracksImpXYZ.push_back(varImpXYZ);
        incjetTracksSignImpXYZ.push_back(varSgnImpXYZ);
        incjetTracksImpXYZSig.push_back(varImpXYZSig);
        incjetTracksSignImpXYZSig.push_back(varSgnImpXYZSig);

        if (jetflavor == JetTaggingSpecies::lightflavour) {
          registry.fill(HIST("h_lfjet_track_pt"), track.pt());
          registry.fill(HIST("h_lfjet_track_phi"), track.phi());
          registry.fill(HIST("h_lfjet_track_eta"), track.eta());
          registry.fill(HIST("h_lfjet_ntracks"), mcdjet.tracks().size());
          registry.fill(HIST("h_lfjet_impact_parameter_xy"), varImpXY);
          registry.fill(HIST("h_lfjet_sign_impact_parameter_xy"), varSgnImpXY);
          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance"), varImpXYSig);
          registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_significance"), varSgnImpXYSig);
          registry.fill(HIST("h_lfjet_pt_sign_impact_parameter_xy_significance"), mcdjet.pt(), varSgnImpXYSig);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz"), varImpXYZ);
          registry.fill(HIST("h_lfjet_sign_impact_parameter_xyz"), varSgnImpXYZ);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance"), varImpXYZSig);
          registry.fill(HIST("h_lfjet_sign_impact_parameter_xyz_significance"), varSgnImpXYZSig);
          lfjetTracksImpXY.push_back(varImpXY);
          lfjetTracksImpXYSig.push_back(varImpXYSig);
          lfjetTracksSignImpXY.push_back(varSgnImpXY);
          lfjetTracksSignImpXYSig.push_back(varSgnImpXYSig);
          lfjetTracksImpXYZ.push_back(varImpXYZ);
          lfjetTracksImpXYZSig.push_back(varImpXYZSig);
        }
        if (jetflavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h_cjet_track_pt"), track.pt());
          registry.fill(HIST("h_cjet_track_phi"), track.phi());
          registry.fill(HIST("h_cjet_track_eta"), track.eta());
          registry.fill(HIST("h_cjet_ntracks"), mcdjet.tracks().size());
          registry.fill(HIST("h_cjet_impact_parameter_xy"), varImpXY);
          registry.fill(HIST("h_cjet_sign_impact_parameter_xy"), varSgnImpXY);
          registry.fill(HIST("h_cjet_impact_parameter_xy_significance"), varImpXYSig);
          registry.fill(HIST("h_cjet_sign_impact_parameter_xy_significance"), varSgnImpXYSig);
          registry.fill(HIST("h_cjet_pt_sign_impact_parameter_xy_significance"), mcdjet.pt(), varSgnImpXYSig);
          registry.fill(HIST("h_cjet_impact_parameter_xyz"), varImpXYZ);
          registry.fill(HIST("h_cjet_sign_impact_parameter_xyz"), varSgnImpXYZ);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_significance"), varImpXYZSig);
          registry.fill(HIST("h_cjet_sign_impact_parameter_xyz_significance"), varSgnImpXYZSig);
          cjetTracksImpXY.push_back(varImpXY);
          cjetTracksImpXYSig.push_back(varImpXYSig);
          cjetTracksSignImpXY.push_back(varSgnImpXY);
          cjetTracksSignImpXYSig.push_back(varSgnImpXYSig);
          cjetTracksImpXYZ.push_back(varImpXYZ);
          cjetTracksImpXYZSig.push_back(varImpXYZSig);
        }
        if (jetflavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h_bjet_track_pt"), track.pt());
          registry.fill(HIST("h_bjet_track_phi"), track.phi());
          registry.fill(HIST("h_bjet_track_eta"), track.eta());
          registry.fill(HIST("h_bjet_ntracks"), mcdjet.tracks().size());
          registry.fill(HIST("h_bjet_impact_parameter_xy"), varImpXY);
          registry.fill(HIST("h_bjet_sign_impact_parameter_xy"), varSgnImpXY);
          registry.fill(HIST("h_bjet_impact_parameter_xy_significance"), varImpXYSig);
          registry.fill(HIST("h_bjet_sign_impact_parameter_xy_significance"), varSgnImpXYSig);
          registry.fill(HIST("h_bjet_pt_sign_impact_parameter_xy_significance"), mcdjet.pt(), varSgnImpXYSig);
          registry.fill(HIST("h_bjet_impact_parameter_xyz"), varImpXYZ);
          registry.fill(HIST("h_bjet_sign_impact_parameter_xyz"), varSgnImpXYZ);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_significance"), varImpXYZSig);
          registry.fill(HIST("h_bjet_sign_impact_parameter_xyz_significance"), varSgnImpXYZSig);
          bjetTracksImpXY.push_back(varImpXY);
          bjetTracksImpXYSig.push_back(varImpXYSig);
          bjetTracksSignImpXY.push_back(varSgnImpXY);
          bjetTracksSignImpXYSig.push_back(varSgnImpXYSig);
          bjetTracksImpXYZ.push_back(varImpXYZ);
          bjetTracksSignImpXYZ.push_back(varSgnImpXYZ);
          bjetTracksImpXYZSig.push_back(varImpXYZSig);
          bjetTracksSignImpXYZSig.push_back(varSgnImpXYZSig);
        }

        // Track counting
        sort(incjetTracksImpXY.begin(), incjetTracksImpXY.end(), std::greater<double>());
        sort(incjetTracksImpXYSig.begin(), incjetTracksImpXYSig.end(), std::greater<double>());
        sort(incjetTracksSignImpXY.begin(), incjetTracksSignImpXY.end(), std::greater<double>());
        sort(incjetTracksSignImpXYSig.begin(), incjetTracksSignImpXYSig.end(), std::greater<double>());
        sort(incjetTracksImpXYZ.begin(), incjetTracksImpXYZ.end(), std::greater<double>());
        sort(incjetTracksImpXYZSig.begin(), incjetTracksImpXYZSig.end(), std::greater<double>());
        if (incjetTracksImpXY.size() > 0) { // N1
          registry.fill(HIST("h_inc_jet_impact_parameter_xy_N1"), incjetTracksImpXY[0]);
          registry.fill(HIST("h_inc_jet_impact_parameter_xy_significance_N1"), incjetTracksImpXYSig[0]);
          registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_N1"), incjetTracksSignImpXY[0]);
          registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_significance_N1"), incjetTracksSignImpXYSig[0]);
          registry.fill(HIST("h_inc_jet_pt_sign_impact_parameter_xy_significance_N1"), mcdjet.pt(), incjetTracksSignImpXYSig[0]);
          registry.fill(HIST("h_inc_jet_impact_parameter_xyz_N1"), incjetTracksImpXYZ[0]);
          registry.fill(HIST("h_inc_jet_impact_parameter_xyz_significance_N1"), incjetTracksImpXYZSig[0]);
        }
        if (incjetTracksImpXY.size() > 1) { // N2
          registry.fill(HIST("h_inc_jet_impact_parameter_xy_N2"), incjetTracksImpXY[1]);
          registry.fill(HIST("h_inc_jet_impact_parameter_xy_significance_N2"), incjetTracksImpXYSig[1]);
          registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_N2"), incjetTracksSignImpXY[1]);
          registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_significance_N2"), incjetTracksSignImpXYSig[1]);
          registry.fill(HIST("h_inc_jet_pt_sign_impact_parameter_xy_significance_N2"), mcdjet.pt(), incjetTracksSignImpXYSig[1]);
          registry.fill(HIST("h_inc_jet_impact_parameter_xyz_N2"), incjetTracksImpXYZ[1]);
          registry.fill(HIST("h_inc_jet_impact_parameter_xyz_significance_N2"), incjetTracksImpXYZSig[1]);
        }
        if (incjetTracksImpXY.size() > 2) { // N3
          registry.fill(HIST("h_inc_jet_impact_parameter_xy_N3"), incjetTracksImpXY[2]);
          registry.fill(HIST("h_inc_jet_impact_parameter_xy_significance_N3"), incjetTracksImpXYSig[2]);
          registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_N3"), incjetTracksSignImpXY[2]);
          registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_significance_N3"), incjetTracksSignImpXYSig[2]);
          registry.fill(HIST("h_inc_jet_pt_sign_impact_parameter_xy_significance_N3"), mcdjet.pt(), incjetTracksSignImpXYSig[2]);
          registry.fill(HIST("h_inc_jet_impact_parameter_xyz_N3"), incjetTracksImpXYZ[2]);
          registry.fill(HIST("h_inc_jet_impact_parameter_xyz_significance_N3"), incjetTracksImpXYZSig[2]);
        }

        if (jetflavor == JetTaggingSpecies::lightflavour && !(lfjetTracksImpXY.empty())) { // lfjet
          sort(lfjetTracksImpXY.begin(), lfjetTracksImpXY.end(), std::greater<double>());
          sort(lfjetTracksImpXYSig.begin(), lfjetTracksImpXYSig.end(), std::greater<double>());
          sort(lfjetTracksSignImpXY.begin(), lfjetTracksSignImpXY.end(), std::greater<double>());
          sort(lfjetTracksSignImpXYSig.begin(), lfjetTracksSignImpXYSig.end(), std::greater<double>());
          sort(lfjetTracksImpXYZ.begin(), lfjetTracksImpXYZ.end(), std::greater<double>());
          sort(lfjetTracksImpXYZSig.begin(), lfjetTracksImpXYZSig.end(), std::greater<double>());
          if (lfjetTracksImpXY.size() > 0) { // N1
            registry.fill(HIST("h_lfjet_impact_parameter_xy_N1"), lfjetTracksImpXY[0]);
            registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_N1"), lfjetTracksImpXYSig[0]);
            registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_N1"), lfjetTracksSignImpXY[0]);
            registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_significance_N1"), lfjetTracksSignImpXYSig[0]);
            registry.fill(HIST("h_lfjet_pt_sign_impact_parameter_xy_significance_N1"), mcdjet.pt(), lfjetTracksSignImpXYSig[0]);
            registry.fill(HIST("h_lfjet_impact_parameter_xyz_N1"), lfjetTracksImpXYZ[0]);
            registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance_N1"), lfjetTracksImpXYZSig[0]);
          }
          if (lfjetTracksImpXY.size() > 1) { // N2
            registry.fill(HIST("h_lfjet_impact_parameter_xy_N2"), lfjetTracksImpXY[1]);
            registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_N2"), lfjetTracksImpXYSig[1]);
            registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_N2"), lfjetTracksSignImpXY[1]);
            registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_significance_N2"), lfjetTracksSignImpXYSig[1]);
            registry.fill(HIST("h_lfjet_pt_sign_impact_parameter_xy_significance_N2"), mcdjet.pt(), lfjetTracksSignImpXYSig[1]);
            registry.fill(HIST("h_lfjet_impact_parameter_xyz_N2"), lfjetTracksImpXYZ[1]);
            registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance_N2"), lfjetTracksImpXYZSig[1]);
          }
          if (lfjetTracksImpXY.size() > 2) { // N3
            registry.fill(HIST("h_lfjet_impact_parameter_xy_N3"), lfjetTracksImpXY[2]);
            registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_N3"), lfjetTracksImpXYSig[2]);
            registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_N3"), lfjetTracksSignImpXY[2]);
            registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_significance_N3"), lfjetTracksSignImpXYSig[2]);
            registry.fill(HIST("h_lfjet_pt_sign_impact_parameter_xy_significance_N3"), mcdjet.pt(), lfjetTracksSignImpXYSig[2]);
            registry.fill(HIST("h_lfjet_impact_parameter_xyz_N3"), lfjetTracksImpXYZ[2]);
            registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance_N3"), lfjetTracksImpXYZSig[2]);
          }
        }
        if (jetflavor == JetTaggingSpecies::charm && !(cjetTracksImpXY.empty())) { // cjet
          sort(cjetTracksImpXY.begin(), cjetTracksImpXY.end(), std::greater<double>());
          sort(cjetTracksImpXYSig.begin(), cjetTracksImpXYSig.end(), std::greater<double>());
          sort(cjetTracksSignImpXY.begin(), cjetTracksSignImpXY.end(), std::greater<double>());
          sort(cjetTracksSignImpXYSig.begin(), cjetTracksSignImpXYSig.end(), std::greater<double>());
          sort(cjetTracksImpXYZ.begin(), cjetTracksImpXYZ.end(), std::greater<double>());
          sort(cjetTracksImpXYZSig.begin(), cjetTracksImpXYZSig.end(), std::greater<double>());
          if (cjetTracksImpXY.size() > 0) { // N1
            registry.fill(HIST("h_cjet_impact_parameter_xy_N1"), cjetTracksImpXY[0]);
            registry.fill(HIST("h_cjet_impact_parameter_xy_significance_N1"), cjetTracksImpXYSig[0]);
            registry.fill(HIST("h_cjet_sign_impact_parameter_xy_N1"), cjetTracksSignImpXY[0]);
            registry.fill(HIST("h_cjet_sign_impact_parameter_xy_significance_N1"), cjetTracksSignImpXYSig[0]);
            registry.fill(HIST("h_cjet_pt_sign_impact_parameter_xy_significance_N1"), mcdjet.pt(), cjetTracksSignImpXYSig[0]);
            registry.fill(HIST("h_cjet_impact_parameter_xyz_N1"), cjetTracksImpXYZ[0]);
            registry.fill(HIST("h_cjet_impact_parameter_xyz_significance_N1"), cjetTracksImpXYZSig[0]);
          }
          if (cjetTracksImpXY.size() > 1) { // N2
            registry.fill(HIST("h_cjet_impact_parameter_xy_N2"), cjetTracksImpXY[1]);
            registry.fill(HIST("h_cjet_impact_parameter_xy_significance_N2"), cjetTracksImpXYSig[1]);
            registry.fill(HIST("h_cjet_sign_impact_parameter_xy_N2"), cjetTracksSignImpXY[1]);
            registry.fill(HIST("h_cjet_sign_impact_parameter_xy_significance_N2"), cjetTracksSignImpXYSig[1]);
            registry.fill(HIST("h_cjet_pt_sign_impact_parameter_xy_significance_N2"), mcdjet.pt(), cjetTracksSignImpXYSig[1]);
            registry.fill(HIST("h_cjet_impact_parameter_xyz_N2"), cjetTracksImpXYZ[1]);
            registry.fill(HIST("h_cjet_impact_parameter_xyz_significance_N2"), cjetTracksImpXYZSig[1]);
          }
          if (cjetTracksImpXY.size() > 2) { // N3
            registry.fill(HIST("h_cjet_impact_parameter_xy_N3"), cjetTracksImpXY[2]);
            registry.fill(HIST("h_cjet_impact_parameter_xy_significance_N3"), cjetTracksImpXYSig[2]);
            registry.fill(HIST("h_cjet_sign_impact_parameter_xy_N3"), cjetTracksSignImpXY[2]);
            registry.fill(HIST("h_cjet_sign_impact_parameter_xy_significance_N3"), cjetTracksSignImpXYSig[2]);
            registry.fill(HIST("h_cjet_pt_sign_impact_parameter_xy_significance_N3"), mcdjet.pt(), cjetTracksSignImpXYSig[2]);
            registry.fill(HIST("h_cjet_impact_parameter_xyz_N3"), cjetTracksImpXYZ[2]);
            registry.fill(HIST("h_cjet_impact_parameter_xyz_significance_N3"), cjetTracksImpXYZSig[2]);
          }
        }
        if (jetflavor == JetTaggingSpecies::beauty && !(bjetTracksImpXY.empty())) { // bjet
          sort(bjetTracksImpXY.begin(), bjetTracksImpXY.end(), std::greater<double>());
          sort(bjetTracksImpXYSig.begin(), bjetTracksImpXYSig.end(), std::greater<double>());
          sort(bjetTracksSignImpXY.begin(), bjetTracksSignImpXY.end(), std::greater<double>());
          sort(bjetTracksSignImpXYSig.begin(), bjetTracksSignImpXYSig.end(), std::greater<double>());
          sort(bjetTracksImpXYZ.begin(), bjetTracksImpXYZ.end(), std::greater<double>());
          sort(bjetTracksImpXYZSig.begin(), bjetTracksImpXYZSig.end(), std::greater<double>());
          if (bjetTracksImpXY.size() > 0) { // N1
            registry.fill(HIST("h_bjet_impact_parameter_xy_N1"), bjetTracksImpXY[0]);
            registry.fill(HIST("h_bjet_impact_parameter_xy_significance_N1"), bjetTracksImpXYSig[0]);
            registry.fill(HIST("h_bjet_sign_impact_parameter_xy_N1"), bjetTracksSignImpXY[0]);
            registry.fill(HIST("h_bjet_sign_impact_parameter_xy_significance_N1"), bjetTracksSignImpXYSig[0]);
            registry.fill(HIST("h_bjet_pt_sign_impact_parameter_xy_significance_N1"), mcdjet.pt(), bjetTracksSignImpXYSig[0]);
            registry.fill(HIST("h_bjet_impact_parameter_xyz_N1"), bjetTracksImpXYZ[0]);
            registry.fill(HIST("h_bjet_impact_parameter_xyz_significance_N1"), bjetTracksImpXYZSig[0]);
          }
          if (bjetTracksImpXY.size() > 1) { // N2
            registry.fill(HIST("h_bjet_impact_parameter_xy_N2"), bjetTracksImpXY[1]);
            registry.fill(HIST("h_bjet_impact_parameter_xy_significance_N2"), bjetTracksImpXYSig[1]);
            registry.fill(HIST("h_bjet_sign_impact_parameter_xy_N2"), bjetTracksSignImpXY[1]);
            registry.fill(HIST("h_bjet_sign_impact_parameter_xy_significance_N2"), bjetTracksSignImpXYSig[1]);
            registry.fill(HIST("h_bjet_pt_sign_impact_parameter_xy_significance_N2"), mcdjet.pt(), bjetTracksSignImpXYSig[1]);
            registry.fill(HIST("h_bjet_impact_parameter_xyz_N2"), bjetTracksImpXYZ[1]);
            registry.fill(HIST("h_bjet_impact_parameter_xyz_significance_N2"), bjetTracksImpXYZSig[1]);
          }
          if (bjetTracksImpXY.size() > 2) { // N3
            registry.fill(HIST("h_bjet_impact_parameter_xy_N3"), bjetTracksImpXY[2]);
            registry.fill(HIST("h_bjet_impact_parameter_xy_significance_N3"), bjetTracksImpXYSig[2]);
            registry.fill(HIST("h_bjet_sign_impact_parameter_xy_N3"), bjetTracksSignImpXY[2]);
            registry.fill(HIST("h_bjet_sign_impact_parameter_xy_significance_N3"), bjetTracksSignImpXYSig[2]);
            registry.fill(HIST("h_bjet_pt_sign_impact_parameter_xy_significance_N3"), mcdjet.pt(), bjetTracksSignImpXYSig[2]);
            registry.fill(HIST("h_bjet_impact_parameter_xyz_N3"), bjetTracksImpXYZ[2]);
            registry.fill(HIST("h_bjet_impact_parameter_xyz_significance_N3"), bjetTracksImpXYZSig[2]);
          }
        }
      }
    }
  }

  template <typename T, typename U, typename V>
  void fillHistogramJPMCD(T const& collision, U const& mcdjets, V const& mcParticles)
  {
    for (const auto& mcdjet : mcdjets) {
      const int jetflavor = mcdjet.origin();
      auto JP = -1.;
      JP = CalculateJP(collision, mcdjet, mcParticles, jetflavor);
      if (JP < 0)
        continue;
      registry.fill(HIST("h_inc_jet_JP"), JP);
      registry.fill(HIST("h_inc_jet_JP_Log"), -1 * TMath::Log(JP));
      registry.fill(HIST("h_inc_jet_pt_JP"), mcdjet.pt(), JP);
      registry.fill(HIST("h_inc_jet_pt_JP_Log"), mcdjet.pt(), -1 * TMath::Log(JP));
      if (jetflavor == JetTaggingSpecies::lightflavour) {
        registry.fill(HIST("h_lfjet_JP"), JP);
        registry.fill(HIST("h_lfjet_JP_Log"), -1 * TMath::Log(JP));
        registry.fill(HIST("h_lfjet_pt_JP"), mcdjet.pt(), JP);
        registry.fill(HIST("h_lfjet_pt_JP_Log"), mcdjet.pt(), -1 * TMath::Log(JP));
      }
      if (jetflavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h_cjet_JP"), JP);
        registry.fill(HIST("h_cjet_JP_Log"), -1 * TMath::Log(JP));
        registry.fill(HIST("h_cjet_pt_JP"), mcdjet.pt(), JP);
        registry.fill(HIST("h_cjet_pt_JP_Log"), mcdjet.pt(), -1 * TMath::Log(JP));
      }
      if (jetflavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_bjet_JP"), JP);
        registry.fill(HIST("h_bjet_JP_Log"), -1 * TMath::Log(JP));
        registry.fill(HIST("h_bjet_pt_JP"), mcdjet.pt(), JP);
        registry.fill(HIST("h_bjet_pt_JP_Log"), mcdjet.pt(), -1 * TMath::Log(JP));
      }
    }
  }

  void processDummy(aod::Collision const& collision)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTC, processDummy, "Dummy process", true);

  void processIPsMCD(aod::Collision const& collision, JetTagTable const& mcdjets, JetTagTracks const& tracks, aod::McParticles const& particles)
  {
    fillHistogramIPsMCD(collision, mcdjets);
  }
  PROCESS_SWITCH(JetTaggerHFTC, processIPsMCD, "Fill impact parameter inpormation for mcd jets", false);

  void processJPMCD(aod::Collision const& collision, JetTagTable const& mcdjets, JetTagTracks const& tracks, aod::McParticles const& particles)
  {
    fillHistogramJPMCD(collision, mcdjets, particles);
  }
  PROCESS_SWITCH(JetTaggerHFTC, processJPMCD, "Fill track and jet probability for mcd jets", false);
};

using JetTaggerTCChargedMCDJets = JetTaggerHFTC<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetTags>>;
// using JetTaggerTCFullMCDJets = JetTaggerHFTC<soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents, aod::FullMCDetectorLevelJetTags>>;
// using JetTaggerTCNeutralMCDJets = JetTaggerHFTC<soa::Join<aod::NeutralMCDetectorLevelJets, aod::NeutralMCDetectorLevelJetConstituents, aod::NeutralMCDetectorLevelJetTags>>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<JetTaggerTCChargedMCDJets>(cfgc,
                                                 SetDefaultProcesses{}, TaskName{"jet-taggerhf-tc-charged"}));
  /*
  tasks.emplace_back(
    adaptAnalysisTask<JetTaggerTCFullMCDJets>(cfgc,
                                         SetDefaultProcesses{}, TaskName{"jet-taggerhf-tc-full"}));

    tasks.emplace_back(
      adaptAnalysisTask<JetTaggerTCNeutralJets>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-taggerhf-tc-neutral"}));
  */
  return WorkflowSpec{tasks};
}
