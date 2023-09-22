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

// \file HfTaggingTask.cxx

// HF trgging task
//
// Authors: hanseo.park@cern.ch

#include "DCAFitter/DCAFitterN.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include <string>
#include "TF1.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;

#include "Framework/runDataProcessing.h"

struct HfTaggingTask {

  static constexpr std::string_view charJetFlavor[] = {"inc_jet", "lfjet", "cjet", "bjet"};

  Configurable<bool> doFitResoFunc{"doFitResoFunc", false, "Do Fit resolution function"};

  // Binning
  ConfigurableAxis binPdgCode{"binPdgCode", {5000, -2500.f, 2500.f}, ""};
  ConfigurableAxis binStatus{"binStatusCode", {200, -99.5f, 100.5f}, ""};
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

  // Axis
  AxisSpec pdgCodeAxis = {binPdgCode, "PDG code"};
  AxisSpec statusAxis = {binStatus, "Status code"};
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

  void init(o2::framework::InitContext&)
  {
    for (int i = 0; i < 4; i++) {
      // common
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
      registry.add(Form("h_%s_particle_pdgcode", charJetFlavor[i].data()), Form("particles' pdg code in %s", charJetFlavor[i].data()), {HistType::kTH1F, {pdgCodeAxis}});
      registry.add(Form("h_%s_particle_statuscode", charJetFlavor[i].data()), Form("particles' status code in %s", charJetFlavor[i].data()), {HistType::kTH1F, {statusAxis}});
      registry.add(Form("h_%s_impact_parameter_xy", charJetFlavor[i].data()), Form("%s impact parameter dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xy", charJetFlavor[i].data()), Form("%s sign impact parameter dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYAxis}});
      registry.add(Form("h_%s_impact_parameter_xy_significance", charJetFlavor[i].data()), Form("%s impact parameter significance dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xy_significance", charJetFlavor[i].data()), Form("%s sign impact parameter significance dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
      registry.add(Form("h_%s_impact_parameter_xyz", charJetFlavor[i].data()), Form("%s impact parameter dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYZAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xyz", charJetFlavor[i].data()), Form("%s sign impact parameter dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYZAxis}});
      registry.add(Form("h_%s_impact_parameter_xyz_significance", charJetFlavor[i].data()), Form("%s impact parameter significance dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xyz_significance", charJetFlavor[i].data()), Form("%s sign impact parameter significance dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});

      for (int j = 1; j < 4; j++) { // Track counting
        registry.add(Form("h_%s_impact_parameter_xy_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xy_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYAxis}});
        registry.add(Form("h_%s_impact_parameter_xy_significance_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter significance dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xy_significance_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
        registry.add(Form("h_%s_impact_parameter_xyz_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYZAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xyz_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYZAxis}});
        registry.add(Form("h_%s_impact_parameter_xyz_significance_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter significance dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xyz_significance_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
      }
      // track and jet probability
      registry.add(Form("h_%s_pos_track_probability", charJetFlavor[i].data()), Form("%s positive sign track probability", charJetFlavor[i].data()), {HistType::kTH1F, {TrackProbabilityAxis}});
      registry.add(Form("h_%s_neg_track_probability", charJetFlavor[i].data()), Form("%s negative sign track probability", charJetFlavor[i].data()), {HistType::kTH1F, {TrackProbabilityAxis}});
      registry.add(Form("h_%s_JP", charJetFlavor[i].data()), Form("%s jet probability", charJetFlavor[i].data()), {HistType::kTH1F, {JetProbabilityAxis}});
      registry.add(Form("h_%s_JP_Log", charJetFlavor[i].data()), Form("Log %s jet probability", charJetFlavor[i].data()), {HistType::kTH1F, {JetProbabilityLogAxis}});
    }

    registry.add("h_jet_entries", "Entries of jets (lf=1, c=2, b=3)", {HistType::kTH1F, {JetEntries}});
    registry.add("h_jet_track_chi2", "track ch2", {HistType::kTH1F, {{100, -10., 10.}}});

    // for Resolution function

    for (int i = 0; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        registry.add(Form("h_%s_sign_impact_parameter_xy_significance_class%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xy} (class%d)", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
      }
    }
  }

  using JetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::StoredTracksIU, aod::TracksCovIU, aod::McTrackLabels>;

  std::unordered_set<int> BHadronSet = {511, 521, 531, 5122, 5112,    // B0, B+, B0s, lambda0b, sigma-b
                                        5212, 5222, 5132, 5232, 5332, // sigma0b, sigma+b, xi-b, xi0b, omeaga-b
                                        5334};                        // omega*-b

  std::unordered_set<int> CHadronSet = {411, 421, 431, 4122, 4112, // D+, D0, D+s, lambda+c, sigma0c
                                        4212, 4222, 4132, 4232, 4332,
                                        4334};

  template <typename T>
  std::pair<int, int> getMotherCodes(T& firstMother)
  {
    std::pair<int, int> motherCodes;
    motherCodes.first = 0;
    motherCodes.second = 0;

    if (!firstMother.has_mothers())
      return motherCodes;
    while (firstMother.has_mothers()) {
      auto mother = firstMother.template mothers_first_as<aod::McParticles>();
      int motherStatusCode = abs(mother.getGenStatusCode());
      if (motherStatusCode == 23 || motherStatusCode == 33 || motherStatusCode == 43 || motherStatusCode == 63) {
        motherCodes.first = motherStatusCode;
        motherCodes.second = abs(mother.pdgCode());
        return motherCodes;
      }
      if (motherStatusCode == 51 && mother.template mothers_first_as<aod::McParticles>().pdgCode() == 21) { // gluon spliting
        while (mother.has_mothers()) {
          auto grandmother = mother.template mothers_first_as<aod::McParticles>();
          int grandMotherStatusCode = abs(grandmother.getGenStatusCode());
          if (grandMotherStatusCode == 23 || grandMotherStatusCode == 33 || grandMotherStatusCode == 43 || grandMotherStatusCode == 63) {
            mother = grandmother;
            break;
          }
          mother = grandmother;
        }
        motherCodes.first = motherStatusCode;
        motherCodes.second = abs(mother.pdgCode());
        return motherCodes;
      }
      firstMother = mother;
    }
    return motherCodes;
  }

  template <typename T>
  bool isBJet(T const& jet)
  {
    std::vector<std::pair<int, int>> candBeauty;
    for (auto& track : jet.template tracks_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::StoredTracksIU, aod::TracksCovIU, aod::McTrackLabels>>()) {
      if (!track.has_mcParticle())
        continue;
      auto const& particle = track.template mcParticle_as<o2::aod::McParticles>();
      if (!particle.has_mothers())
        continue;
      auto mother = particle;
      if (!mother.has_mothers())
        continue;
      while (mother.has_mothers()) {
        auto pdgcode = abs(mother.pdgCode());
        if (BHadronSet.count(pdgcode)) { // bjet
          std::pair<int, int> motherCodes = getMotherCodes(mother);
          candBeauty.push_back(std::make_pair(motherCodes.first, motherCodes.second));
        }
        auto grandmother = mother.template mothers_first_as<aod::McParticles>();
        mother = grandmother;
      }
    }
    for (const auto& cand : candBeauty) { // remove gluon spliting effect
      if (!(cand.first == 51 && cand.second == 21))
        return true;
    }
    return false;
  }

  template <typename T>
  bool isCJet(T const& jet)
  {
    std::vector<std::pair<int, int>> candCharm;
    for (auto& track : jet.template tracks_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::StoredTracksIU, aod::TracksCovIU, aod::McTrackLabels>>()) {
      if (!track.has_mcParticle())
        continue;
      auto const& particle = track.template mcParticle_as<o2::aod::McParticles>();
      auto mother = particle;
      if (!mother.has_mothers())
        continue;
      while (mother.has_mothers()) {
        auto pdgcode = abs(mother.pdgCode());
        if (CHadronSet.count(pdgcode)) { // cjet
          std::pair<int, int> motherCodes = getMotherCodes(mother);
          candCharm.push_back(std::make_pair(motherCodes.first, motherCodes.second));
        }
        auto grandmother = mother.template mothers_first_as<aod::McParticles>();
        mother = grandmother;
      }
    }
    for (const auto& cand : candCharm) { // remove gluon spliting effect
      if (!(cand.first == 51 && cand.second == 21))
        return true;
    }
    return false;
  }

  template <typename T>
  int getJetFlavor(T const& jet)
  {
    bool isB = isBJet(jet);
    bool isC = isCJet(jet);
    if (isB)
      return 3; // b-jet
    if (!isB && isC)
      return 2; // c-jet
    if (!(isB || isC))
      return 1; // lf-jet
    return 0;
  }
  template <typename T>
  bool trackSelectionJP(T const& jet)
  {
    if (jet.tracks().size() < 2)
      return 0;

    return 1;
  }

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

  template <typename T, typename U, typename V, typename D = double>
  double CalculateJP(T const& collision, U const& jet, V const& mcParticles, D& jetflavor)
  {
    double JP = -1.;
    if (!trackSelectionJP(jet))
      return -1.;
    SetFitResoFunc();
    std::vector<double> jetTracksPt;

    double trackjetProb = 1.;
    if (jetflavor == 0)
      return -1;
    for (auto& track : jet.template tracks_as<JetTracks>()) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for track, skip...");
        continue;
      }
      auto particle = track.mcParticle();
      if (!particle.has_mothers()) {
        LOGF(warning, "No mother particle for particle, skip...");
        continue;
      }

      double sign, varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig, varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig;
      SetSignImpactParameterSignificance(collision, jet, track, sign, varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig, varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig);

      if (jetflavor == 0)
        continue;

      double ProbTrack = 0.;

      auto minSignImpXYSig = -10;
      ProbTrack = fResoFuncincjet->Integral(minSignImpXYSig, -1 * TMath::Abs(varSignImpXYSig)) / fResoFuncincjet->Integral(minSignImpXYSig, 0);
      if (TMath::Abs(varSignImpXYSig) > TMath::Abs(minSignImpXYSig))
        varSignImpXYSig = 9.99; // Limit to function definition range
      if (sign > 0)
        registry.fill(HIST("h_inc_jet_pos_track_probability"), ProbTrack);
      else
        registry.fill(HIST("h_inc_jet_neg_track_probability"), ProbTrack);
      if (jetflavor == 1) { // lf
        ProbTrack = fResoFunclfjet->Integral(minSignImpXYSig, -1 * TMath::Abs(varSignImpXYSig)) / fResoFunclfjet->Integral(minSignImpXYSig, 0);
        if (sign > 0)
          registry.fill(HIST("h_lfjet_pos_track_probability"), ProbTrack);
        else
          registry.fill(HIST("h_lfjet_neg_track_probability"), ProbTrack);
      }
      if (jetflavor == 2) { // charm
        ProbTrack = fResoFunccjet->Integral(minSignImpXYSig, -1 * TMath::Abs(varSignImpXYSig)) / fResoFunccjet->Integral(minSignImpXYSig, 0);
        if (sign > 0)
          registry.fill(HIST("h_cjet_pos_track_probability"), ProbTrack);
        else
          registry.fill(HIST("h_cjet_neg_track_probability"), ProbTrack);
      }
      if (jetflavor == 3) { // beauty
        ProbTrack = fResoFuncbjet->Integral(minSignImpXYSig, -1 * TMath::Abs(varSignImpXYSig)) / fResoFuncbjet->Integral(minSignImpXYSig, 0);
        if (sign > 0)
          registry.fill(HIST("h_bjet_pos_track_probability"), ProbTrack);
        else
          registry.fill(HIST("h_bjet_neg_track_probability"), ProbTrack);
      }
      if (sign > 0) { // only take postive sign track
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

  template <typename T, typename U, typename V>
  void SetSignImpactParameterSignificance(T const& collision, U const& jet, V const& track, double& sign, double& IPxy, double& SignIPxy, double& IPxySig, double& SignIPxySig, double& IPxyz, double& SignIPxyz, double& IPxyzSig, double& SignIPxyzSig)
  {
    sign = TMath::Sign(1, (track.x() - collision.posX()) * jet.px() + (track.y() - collision.posY()) * jet.py() + (track.z() - collision.posZ()) * jet.pz());
    IPxy = track.dcaXY();
    SignIPxy = sign * TMath::Abs(IPxy);
    IPxySig = IPxy / TMath::Sqrt(track.sigmaDcaXY2());
    SignIPxySig = sign * TMath::Abs(IPxySig);

    IPxyz = TMath::Sqrt(track.dcaXY() * track.dcaXY() + track.dcaZ() * track.dcaZ()); // TODO: track.getSigmaZY()
    SignIPxyz = sign * TMath::Abs(IPxyz);
    IPxyzSig = IPxyz / TMath::Sqrt(TMath::Sqrt(track.sigmaDcaXY2() * track.sigmaDcaXY2() + track.sigmaDcaZ2() * track.sigmaDcaZ2())); // TODO: cal cov
    SignIPxyzSig = sign * TMath::Abs(IPxyzSig);
  }

  template <typename T>
  int getQualityTrack(T const& track)
  { // TODO

    if (!track.hasITS())
      return 0;
    if (track.itsChi2NCl() > 2)
      return 1;

    if (track.itsChi2NCl() < 2 && track.pt() < 2)
      return 2;
    if (track.itsChi2NCl() < 2 && track.pt() > 2)
      return 3;

    return 0;
  }

  template <typename T, typename U, typename V>
  void fillIPsMCDHistograms(T const& collision, U const& jets, V const& mcParticles, float weight = 1.0)
  {
    for (const auto& jet : jets) {
      // common var
      std::vector<double> incjetTracksImpXY, lfjetTracksImpXY, cjetTracksImpXY, bjetTracksImpXY;
      std::vector<double> incjetTracksSignImpXY, lfjetTracksSignImpXY, cjetTracksSignImpXY, bjetTracksSignImpXY;
      std::vector<double> incjetTracksImpXYSig, lfjetTracksImpXYSig, cjetTracksImpXYSig, bjetTracksImpXYSig;
      std::vector<double> incjetTracksSignImpXYSig, lfjetTracksSignImpXYSig, cjetTracksSignImpXYSig, bjetTracksSignImpXYSig;
      std::vector<double> incjetTracksImpXYZ, lfjetTracksImpXYZ, cjetTracksImpXYZ, bjetTracksImpXYZ;
      std::vector<double> incjetTracksSignImpXYZ, lfjetTracksSignImpXYZ, cjetTracksSignImpXYZ, bjetTracksSignImpXYZ;
      std::vector<double> incjetTracksImpXYZSig, lfjetTracksImpXYZSig, cjetTracksImpXYZSig, bjetTracksImpXYZSig;
      std::vector<double> incjetTracksSignImpXYZSig, lfjetTracksSignImpXYZSig, cjetTracksSignImpXYZSig, bjetTracksSignImpXYZSig;

      const int jetflavor = getJetFlavor(jet);
      registry.fill(HIST("h_inc_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_inc_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_inc_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_inc_jet_ntracks"), jet.tracks().size(), weight);

      if (jetflavor == 1) { // lfjet
        registry.fill(HIST("h_lfjet_pt"), jet.pt(), weight);
        registry.fill(HIST("h_lfjet_eta"), jet.eta(), weight);
        registry.fill(HIST("h_lfjet_phi"), jet.phi(), weight);
        registry.fill(HIST("h_lfjet_ntracks"), jet.tracks().size(), weight);
      }
      if (jetflavor == 2) { // cjet
        registry.fill(HIST("h_cjet_pt"), jet.pt(), weight);
        registry.fill(HIST("h_cjet_eta"), jet.eta(), weight);
        registry.fill(HIST("h_cjet_phi"), jet.phi(), weight);
        registry.fill(HIST("h_cjet_ntracks"), jet.tracks().size(), weight);
      }
      if (jetflavor == 3) { // bjet
        registry.fill(HIST("h_bjet_pt"), jet.pt(), weight);
        registry.fill(HIST("h_bjet_eta"), jet.eta(), weight);
        registry.fill(HIST("h_bjet_phi"), jet.phi(), weight);
        registry.fill(HIST("h_bjet_ntracks"), jet.tracks().size(), weight);
      }
      for (const auto& track : jet.template tracks_as<JetTracks>()) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC particle for track, skip...");
          continue;
        }

        if (!track.isPrimaryTrack()) {
          LOGF(warning, "No primary track, skip...");
          continue;
        }

        double sign, varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig, varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig;
        SetSignImpactParameterSignificance(collision, jet, track, sign, varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig, varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig);
        registry.fill(HIST("h_inc_jet_track_pt"), track.pt());
        registry.fill(HIST("h_inc_jet_track_phi"), track.phi());
        registry.fill(HIST("h_inc_jet_track_eta"), track.eta());
        registry.fill(HIST("h_inc_jet_impact_parameter_xy"), varImpXY);
        registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy"), varSignImpXY);
        registry.fill(HIST("h_inc_jet_impact_parameter_xy_significance"), varImpXYSig);
        registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_significance"), varSignImpXYSig);
        registry.fill(HIST("h_inc_jet_impact_parameter_xyz"), varImpXYZ);
        registry.fill(HIST("h_inc_jet_impact_parameter_xyz_significance"), varImpXYZSig);
        registry.fill(HIST("h_inc_jet_particle_pdgcode"), track.mcParticle().pdgCode());
        registry.fill(HIST("h_inc_jet_particle_statuscode"), track.mcParticle().getGenStatusCode());
        incjetTracksImpXY.push_back(varImpXY);
        incjetTracksImpXYSig.push_back(varImpXYSig);
        incjetTracksSignImpXY.push_back(varSignImpXY);
        incjetTracksSignImpXYSig.push_back(varSignImpXYSig);
        incjetTracksImpXYZ.push_back(varImpXYZ);
        incjetTracksImpXYZSig.push_back(varImpXYZSig);
        if (jetflavor == 1) { // lf
          registry.fill(HIST("h_lfjet_track_pt"), track.pt());
          registry.fill(HIST("h_lfjet_track_phi"), track.phi());
          registry.fill(HIST("h_lfjet_track_eta"), track.eta());
          registry.fill(HIST("h_lfjet_ntracks"), jet.tracks().size(), weight);
          registry.fill(HIST("h_lfjet_impact_parameter_xy"), varImpXY);
          registry.fill(HIST("h_lfjet_sign_impact_parameter_xy"), varSignImpXY);
          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance"), varImpXYSig);
          registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_significance"), varSignImpXYSig);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz"), varImpXYZ);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance"), varImpXYZSig);
          registry.fill(HIST("h_lfjet_particle_pdgcode"), track.mcParticle().pdgCode());
          registry.fill(HIST("h_lfjet_particle_statuscode"), track.mcParticle().getGenStatusCode());
          lfjetTracksImpXY.push_back(varImpXY);
          lfjetTracksImpXYSig.push_back(varImpXYSig);
          lfjetTracksSignImpXY.push_back(varSignImpXY);
          lfjetTracksSignImpXYSig.push_back(varSignImpXYSig);
          lfjetTracksImpXYZ.push_back(varImpXYZ);
          lfjetTracksImpXYZSig.push_back(varImpXYZSig);
        }
        if (jetflavor == 2) { // charm
          registry.fill(HIST("h_cjet_track_pt"), track.pt());
          registry.fill(HIST("h_cjet_track_phi"), track.phi());
          registry.fill(HIST("h_cjet_track_eta"), track.eta());
          registry.fill(HIST("h_cjet_ntracks"), jet.tracks().size(), weight);
          registry.fill(HIST("h_cjet_impact_parameter_xy"), varImpXY);
          registry.fill(HIST("h_cjet_sign_impact_parameter_xy"), varSignImpXY);
          registry.fill(HIST("h_cjet_impact_parameter_xy_significance"), varImpXYSig);
          registry.fill(HIST("h_cjet_sign_impact_parameter_xy_significance"), varSignImpXYSig);
          registry.fill(HIST("h_cjet_impact_parameter_xyz"), varImpXYZ);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_significance"), varImpXYZSig);
          cjetTracksImpXY.push_back(varImpXY);
          cjetTracksImpXYSig.push_back(varImpXYSig);
          cjetTracksSignImpXY.push_back(varSignImpXY);
          cjetTracksSignImpXYSig.push_back(varSignImpXYSig);
          cjetTracksImpXYZ.push_back(varImpXYZ);
          cjetTracksImpXYZSig.push_back(varImpXYZSig);
        }
        if (jetflavor == 3) { // beauty
          registry.fill(HIST("h_bjet_track_pt"), track.pt());
          registry.fill(HIST("h_bjet_track_phi"), track.phi());
          registry.fill(HIST("h_bjet_track_eta"), track.eta());
          registry.fill(HIST("h_bjet_ntracks"), jet.tracks().size(), weight);
          registry.fill(HIST("h_bjet_impact_parameter_xy"), varImpXY);
          registry.fill(HIST("h_bjet_sign_impact_parameter_xy"), varSignImpXY);
          registry.fill(HIST("h_bjet_impact_parameter_xy_significance"), varImpXYSig);
          registry.fill(HIST("h_bjet_sign_impact_parameter_xy_significance"), varSignImpXYSig);
          registry.fill(HIST("h_bjet_impact_parameter_xyz"), varImpXYZ);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_significance"), varImpXYZSig);
          bjetTracksImpXY.push_back(varImpXY);
          bjetTracksImpXYSig.push_back(varImpXYSig);
          bjetTracksSignImpXY.push_back(varSignImpXY);
          bjetTracksSignImpXYSig.push_back(varSignImpXYSig);
          bjetTracksImpXYZ.push_back(varImpXYZ);
          bjetTracksImpXYZSig.push_back(varImpXYZSig);
        }
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
        registry.fill(HIST("h_inc_jet_impact_parameter_xyz_N1"), incjetTracksImpXYZ[0]);
        registry.fill(HIST("h_inc_jet_impact_parameter_xyz_significance_N1"), incjetTracksImpXYZSig[0]);
      }
      if (incjetTracksImpXY.size() > 1) { // N2
        registry.fill(HIST("h_inc_jet_impact_parameter_xy_N2"), incjetTracksImpXY[1]);
        registry.fill(HIST("h_inc_jet_impact_parameter_xy_significance_N2"), incjetTracksImpXYSig[1]);
        registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_N2"), incjetTracksSignImpXY[1]);
        registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_significance_N2"), incjetTracksSignImpXYSig[1]);
        registry.fill(HIST("h_inc_jet_impact_parameter_xyz_N2"), incjetTracksImpXYZ[1]);
        registry.fill(HIST("h_inc_jet_impact_parameter_xyz_significance_N2"), incjetTracksImpXYZSig[1]);
      }
      if (incjetTracksImpXY.size() > 2) { // N3
        registry.fill(HIST("h_inc_jet_impact_parameter_xy_N3"), incjetTracksImpXY[2]);
        registry.fill(HIST("h_inc_jet_impact_parameter_xy_significance_N3"), incjetTracksImpXYSig[2]);
        registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_N3"), incjetTracksSignImpXY[2]);
        registry.fill(HIST("h_inc_jet_sign_impact_parameter_xy_significance_N3"), incjetTracksSignImpXYSig[2]);
        registry.fill(HIST("h_inc_jet_impact_parameter_xyz_N3"), incjetTracksImpXYZ[2]);
        registry.fill(HIST("h_inc_jet_impact_parameter_xyz_significance_N3"), incjetTracksImpXYZSig[2]);
      }

      if (jetflavor == 1 && !(lfjetTracksImpXY.empty())) { // lfjet
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
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_N1"), lfjetTracksImpXYZ[0]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance_N1"), lfjetTracksImpXYZSig[0]);
        }
        if (lfjetTracksImpXY.size() > 1) { // N2
          registry.fill(HIST("h_lfjet_impact_parameter_xy_N2"), lfjetTracksImpXY[1]);
          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_N2"), lfjetTracksImpXYSig[1]);
          registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_N2"), lfjetTracksSignImpXY[1]);
          registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_significance_N2"), lfjetTracksSignImpXYSig[1]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_N2"), lfjetTracksImpXYZ[1]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance_N2"), lfjetTracksImpXYZSig[1]);
        }
        if (lfjetTracksImpXY.size() > 2) { // N3
          registry.fill(HIST("h_lfjet_impact_parameter_xy_N3"), lfjetTracksImpXY[2]);
          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_N3"), lfjetTracksImpXYSig[2]);
          registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_N3"), lfjetTracksSignImpXY[2]);
          registry.fill(HIST("h_lfjet_sign_impact_parameter_xy_significance_N3"), lfjetTracksSignImpXYSig[2]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_N3"), lfjetTracksImpXYZ[2]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance_N3"), lfjetTracksImpXYZSig[2]);
        }
      }
      if (jetflavor == 2 && !(cjetTracksImpXY.empty())) { // cjet
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
          registry.fill(HIST("h_cjet_impact_parameter_xyz_N1"), cjetTracksImpXYZ[0]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_significance_N1"), cjetTracksImpXYZSig[0]);
        }
        if (cjetTracksImpXY.size() > 1) { // N2
          registry.fill(HIST("h_cjet_impact_parameter_xy_N2"), cjetTracksImpXY[1]);
          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_N2"), cjetTracksImpXYSig[1]);
          registry.fill(HIST("h_cjet_sign_impact_parameter_xy_N2"), cjetTracksSignImpXY[1]);
          registry.fill(HIST("h_cjet_sign_impact_parameter_xy_significance_N2"), cjetTracksSignImpXYSig[1]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_N2"), cjetTracksImpXYZ[1]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_significance_N2"), cjetTracksImpXYZSig[1]);
        }
        if (cjetTracksImpXY.size() > 2) { // N3
          registry.fill(HIST("h_cjet_impact_parameter_xy_N3"), cjetTracksImpXY[2]);
          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_N3"), cjetTracksImpXYSig[2]);
          registry.fill(HIST("h_cjet_sign_impact_parameter_xy_N3"), cjetTracksSignImpXY[2]);
          registry.fill(HIST("h_cjet_sign_impact_parameter_xy_significance_N3"), cjetTracksSignImpXYSig[2]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_N3"), cjetTracksImpXYZ[2]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_significance_N3"), cjetTracksImpXYZSig[2]);
        }
      }
      if (jetflavor == 3 && !(bjetTracksImpXY.empty())) { // bjet
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
          registry.fill(HIST("h_bjet_impact_parameter_xyz_N1"), bjetTracksImpXYZ[0]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_significance_N1"), bjetTracksImpXYZSig[0]);
        }
        if (bjetTracksImpXY.size() > 1) { // N2
          registry.fill(HIST("h_bjet_impact_parameter_xy_N2"), bjetTracksImpXY[1]);
          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_N2"), bjetTracksImpXYSig[1]);
          registry.fill(HIST("h_bjet_sign_impact_parameter_xy_N2"), bjetTracksSignImpXY[1]);
          registry.fill(HIST("h_bjet_sign_impact_parameter_xy_significance_N2"), bjetTracksSignImpXYSig[1]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_N2"), bjetTracksImpXYZ[1]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_significance_N2"), bjetTracksImpXYZSig[1]);
        }
        if (bjetTracksImpXY.size() > 2) { // N3
          registry.fill(HIST("h_bjet_impact_parameter_xy_N3"), bjetTracksImpXY[2]);
          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_N3"), bjetTracksImpXYSig[2]);
          registry.fill(HIST("h_bjet_sign_impact_parameter_xy_N3"), bjetTracksSignImpXY[2]);
          registry.fill(HIST("h_bjet_sign_impact_parameter_xy_significance_N3"), bjetTracksSignImpXYSig[2]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_N3"), bjetTracksImpXYZ[2]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_significance_N3"), bjetTracksImpXYZSig[2]);
        }
      }
    }
  }

  template <typename T, typename U, typename V>
  void fillJetProbMCDHistograms(T const& collision, U const& jets, V const& mcParticles, float weight = 1.0)
  {
    for (const auto& jet : jets) {

      const int jetflavor = getJetFlavor(jet);
      double JP = -1.;
      JP = CalculateJP(collision, jet, mcParticles, jetflavor);
      if (JP < 0)
        continue;
      if (jetflavor == 1) { // lf
        registry.fill(HIST("h_lfjet_JP"), JP);
        registry.fill(HIST("h_lfjet_JP_Log"), -1 * TMath::Log(JP));
      }
      if (jetflavor == 2) { // charm
        registry.fill(HIST("h_cjet_JP"), JP);
        registry.fill(HIST("h_cjet_JP_Log"), -1 * TMath::Log(JP));
      }
      if (jetflavor == 3) { // beauty
        registry.fill(HIST("h_bjet_JP"), JP);
        registry.fill(HIST("h_bjet_JP_Log"), -1 * TMath::Log(JP));
      }
    }
  }

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(HfTaggingTask, processDummy, "Dummy process function turned on by default", true);

  void processIPsJetsMCD(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                         soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                         JetTracks const& tracks,
                         aod::McCollisions const& mcCollisions,
                         aod::McParticles const& mcParticles)
  {
    fillIPsMCDHistograms(collision, jets, mcParticles);
  }
  PROCESS_SWITCH(HfTaggingTask, processIPsJetsMCD, "Task of impact parameter for heavy flavor tagging flavor (MCD)", false);

  void processJetProbMCD(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                         soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                         JetTracks const& tracks,
                         aod::McCollisions const& mcCollisions,
                         aod::McParticles const& mcParticles)
  {
    fillJetProbMCDHistograms(collision, jets, mcParticles);
  }
  PROCESS_SWITCH(HfTaggingTask, processJetProbMCD, "Task of jet probability for heavy flavor tagging (MCD)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<HfTaggingTask>(cfgc, TaskName{"hf-tagging-task"})}; }
