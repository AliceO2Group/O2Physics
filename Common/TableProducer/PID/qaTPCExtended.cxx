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
/// \file qaTPCExtended.cxx
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>, GSI Darmstadt
/// \brief Task to give QA output for TPC PID response based on external classes (TOF cut, V0s, etc)

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/StrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using TPCV0Tracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi,
                              aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TrackSelection>;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{         //runtime customisation goes here
  {"useV0", VariantType::Int, 0, {"Use V0 information for QA"}}
  };
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

using namespace o2::dataformats;

struct QaTpcTof {
  //Configurables
  Configurable<float> cutTOF{"cutTOF", 3.f, "TOF nsigma cut for TPC-TOF PID"};
  Configurable<int> pBins{"pBins", 400, "Number of momentum bins"};
  Configurable<float> pMin{"pMin", 0.01, "Lower limit in momentum"};
  Configurable<float> pMax{"pMax", 20., "Upper limit in momentum"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 200, "Number of bins for TPC nSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.f, "Lower limit for TPC nSigma"};
  Configurable<float> maxNSigma{"maxNSigma", 10.f, "Upper limit for TPC nSigma"};

  HistogramRegistry hists{"HistogramsTPCPIDQA", {}, OutputObjHandlingPolicy::QAObject};

  static constexpr int Np = 9;
  static constexpr std::string_view hnsigmaTPC[Np] = {"nsigmaTPC/El", "nsigmaTPC/Mu", "nsigmaTPC/Pi",
                                                      "nsigmaTPC/Ka", "nsigmaTPC/Pr", "nsigmaTPC/De",
                                                      "nsigmaTPC/Tr", "nsigmaTPC/He", "nsigmaTPC/Al"};
  static constexpr std::string_view hnsigmaTPCTOF[Np] = {"nsigmaTPCTOF/El", "nsigmaTPCTOF/Mu", "nsigmaTPCTOF/Pi",
                                                         "nsigmaTPCTOF/Ka", "nsigmaTPCTOF/Pr", "nsigmaTPCTOF/De",
                                                         "nsigmaTPCTOF/Tr", "nsigmaTPCTOF/He", "nsigmaTPCTOF/Al"};

  static constexpr std::string_view hnsigmaTOF[Np] = {"nsigmaTOF/El", "nsigmaTOF/Mu", "nsigmaTOF/Pi",
                                                      "nsigmaTOF/Ka", "nsigmaTOF/Pr", "nsigmaTOF/De",
                                                      "nsigmaTOF/Tr", "nsigmaTOF/He", "nsigmaTOF/Al"};
  static constexpr std::string_view hnsigmaTOFAfter[Np] = {"nsigmaTOFAfter/El", "nsigmaTOFAfter/Mu", "nsigmaTOFAfter/Pi",
                                                           "nsigmaTOFAfter/Ka", "nsigmaTOFAfter/Pr", "nsigmaTOFAfter/De",
                                                           "nsigmaTOFAfter/Tr", "nsigmaTOFAfter/He", "nsigmaTOFAfter/Al"};
  static constexpr const char* partName[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};

  template <uint8_t i>
  void addTPCQAParticleHistos()
  {
    AxisSpec pAxis{pBins, pMin, pMax, "#it{p} [GeV/#it{c}]"};
    pAxis.makeLogaritmic();
    AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, "TPC n_{#sigma}"};
    AxisSpec tofnSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, "TOF n_{#sigma}"};

    hists.add(hnsigmaTPC[i].data(), Form("TPC signal (%s) without TOF cut", partName[i]), kTH2F, {pAxis, nSigmaAxis});
    hists.add(hnsigmaTPCTOF[i].data(), Form("TPC signal (%s) after %.2f#sigma TOF cut", partName[i], double(cutTOF)), kTH2F, {pAxis, nSigmaAxis});
    hists.add(hnsigmaTOF[i].data(), "TOF signal", kTH2F, {pAxis, tofnSigmaAxis});
    hists.add(hnsigmaTOFAfter[i].data(), "TOF signal after TOF cut", kTH2F, {pAxis, tofnSigmaAxis});

  } //addParticleHistos

  void init(InitContext&)
  {

    addTPCQAParticleHistos<0>();
    addTPCQAParticleHistos<1>();
    addTPCQAParticleHistos<2>();
    addTPCQAParticleHistos<3>();
    addTPCQAParticleHistos<4>();
    addTPCQAParticleHistos<5>();
    addTPCQAParticleHistos<6>();
    addTPCQAParticleHistos<7>();
    addTPCQAParticleHistos<8>();

  } //init

  template <uint8_t i, typename T>
  void fillTPCQAParticleHistos(const T& t, const float mom, const float tofNSigma, const float tpcNSigma)
  {
    // Fill TPC-TOF histograms before/after nsigma cut on TOF
    hists.fill(HIST(hnsigmaTPC[i]), mom, tpcNSigma);
    if (abs(tofNSigma) < cutTOF)
      hists.fill(HIST(hnsigmaTPCTOF[i]), mom, tpcNSigma);

    hists.fill(HIST(hnsigmaTOF[i]), mom, tofNSigma);
    if (abs(tofNSigma) < cutTOF)
      hists.fill(HIST(hnsigmaTOFAfter[i]), mom, tofNSigma);
  } //fillParticleHistos

  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra,
                                                          aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                                          aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                                          aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl,
                                                          aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                                          aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe,
                                                          aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                                                          aod::TrackSelection> const& tracks)
  {
    for (auto t : tracks) {
      const float mom = t.tpcInnerParam();
      fillTPCQAParticleHistos<0>(t, mom, t.tofNSigmaEl(), t.tpcNSigmaEl());
      fillTPCQAParticleHistos<1>(t, mom, t.tofNSigmaMu(), t.tpcNSigmaMu());
      fillTPCQAParticleHistos<2>(t, mom, t.tofNSigmaPi(), t.tpcNSigmaPi());
      fillTPCQAParticleHistos<3>(t, mom, t.tofNSigmaKa(), t.tpcNSigmaKa());
      fillTPCQAParticleHistos<4>(t, mom, t.tofNSigmaPr(), t.tpcNSigmaPr());
      fillTPCQAParticleHistos<5>(t, mom, t.tofNSigmaDe(), t.tpcNSigmaDe());
      fillTPCQAParticleHistos<6>(t, mom, t.tofNSigmaTr(), t.tpcNSigmaTr());
      fillTPCQAParticleHistos<7>(t, mom, t.tofNSigmaHe(), t.tpcNSigmaHe());
      fillTPCQAParticleHistos<8>(t, mom, t.tofNSigmaAl(), t.tpcNSigmaAl());

    } //for
  }   //process
};    //struct QaTPCPID

struct QaTpcV0 {
  static constexpr int NpV0 = 4;
  enum EV0DaughID { kEle = 0,
                    kPi = 1,
                    kKa = 2,
                    kPr = 3 };
  static constexpr const char* partName[NpV0] = {"e", "#pi", "K", "p"};
  static constexpr std::string_view hdEdxV0[NpV0] = {"dEdxV0/El", "dEdxV0/Pi",
                                                     "dEdxV0/Ka", "dEdxV0/Pr"};
  static constexpr std::string_view hdEdxDiffV0[NpV0] = {"dEdxDiffTPCV0/El", "dEdxDiffTPCV0/Pi",
                                                         "dEdxDiffTPCV0/Ka", "dEdxDiffTPCV0/Pr"};
  static constexpr std::string_view hnsigmaV0[NpV0] = {"nsigmaTPCV0/El", "nsigmaTPCV0/Pi",
                                                       "nsigmaTPCV0/Ka", "nsigmaTPCV0/Pr"};
  static constexpr std::string_view hnsigmaV0VsEta[NpV0] = {"nsigmaTPCV0VsEta/El", "nsigmaTPCV0VsEta/Pi",
                                                       "nsigmaTPCV0VsEta/Ka", "nsigmaTPCV0VsEta/Pr"};
  HistogramRegistry histos{"TPCPIDQA_V0", {}, OutputObjHandlingPolicy::QAObject};
  Configurable<int> nBinsP{"nBinsP", 400, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 200, "Number of bins for TPC nSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.f, "Lower limit for TPC nSigma"};
  Configurable<float> maxNSigma{"maxNSigma", 10.f, "Upper limit for TPC nSigma"};

  //Definition of V0 preselection cuts for K0S, Lambda, Anti-lambda
  //K0S
  static constexpr const float cutQTK0[2] = {0.1075f, 0.215f};
  static constexpr const float cutAPK0[2] = {0.199f, 0.8f};
  //Lambda/Antilambda
  static constexpr const float cutQTL = 0.03f;
  static constexpr const float cutAlphaL[2] = {0.35f, 0.7f};
  static constexpr const float cutAlphaAL[2] = {-0.7f, -0.35f};
  static constexpr const float cutAPL[3] = {0.107f, -0.69f, 0.5f};

  template <uint8_t i>
  void addV0Histos()
  {
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    pAxis.makeLogaritmic();
    

    // corrected dE/dx for clean V0
    AxisSpec expAxis{1000, 0, 1000, Form("d#it{E}/d#it{x}_(%s from V^{0}) A.U.", partName[i])};
    histos.add(hdEdxV0[i].data(), "", kTH2F, {pAxis, expAxis});

    AxisSpec deltaAxis{1000, -500, 500, Form("d#it{E}/d#it{x} - d#it{E}/d#it{x}(expected(%s from V^{0})", partName[i])};
    histos.add(hdEdxDiffV0[i].data(), "", kTH2F, {pAxis, deltaAxis});

    AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, Form("n_{#sigma}^{TPC}(%s from V^{0})", partName[i])};
    histos.add(hnsigmaV0[i].data(), "", kTH2F, {pAxis, nSigmaAxis});

    AxisSpec etaAxis{1000, -1, 1, "#eta"};
    histos.add(hnsigmaV0VsEta[i].data(), "", kTH2F, {etaAxis, nSigmaAxis});

  } //addV0Histos

  template <uint8_t i, typename T>
  void fillV0Histos(const T& t, const float mom, const float exp_diff, const float nsigma)
  {
    histos.fill(HIST(hdEdxV0[i]), mom, t.tpcSignal() - exp_diff);
    histos.fill(HIST(hdEdxDiffV0[i]), mom, exp_diff);
    histos.fill(HIST(hnsigmaV0[i]), t.p(), nsigma);
    histos.fill(HIST(hnsigmaV0VsEta[i]),t.eta(),nsigma);
  }

  void init(o2::framework::InitContext&)
  {
    addV0Histos<0>();
    addV0Histos<1>();
    addV0Histos<2>();
    addV0Histos<3>();

  } //init

  void process(aod::Collision const& collision, aod::V0Datas const& v0s, TPCV0Tracks const& tracks)
  {
    for (auto v0 : v0s) { //for loop on built v0 candidates
      // initialise dynamic variables
      float alpha = v0.alpha();
      float qt = v0.qtarm();
      auto posTrack = v0.posTrack_as<TPCV0Tracks>();
      auto negTrack = v0.negTrack_as<TPCV0Tracks>();
      // Check for K0
      float q = cutAPK0[0] * sqrt(abs(1 - alpha * alpha / (cutAPK0[1] * cutAPK0[1])));
      if ((qt > cutQTK0[0]) && (qt < cutQTK0[1]) && (qt > q)) {
        // Treat as K0 (both tracks pions)
        fillV0Histos<kPi>(posTrack, posTrack.tpcInnerParam(), posTrack.tpcExpSignalDiffPi(), posTrack.tpcNSigmaPi());
        fillV0Histos<kPi>(negTrack, negTrack.tpcInnerParam(), negTrack.tpcExpSignalDiffPi(), negTrack.tpcNSigmaPi());
      }

      // Check for Lambda
      q = cutAPL[0] * sqrt(abs(1 - ((alpha + cutAPL[1]) * (alpha + cutAPL[1])) / (cutAPL[2] * cutAPL[2])));
      if ((alpha > cutAlphaL[0]) && (alpha < cutAlphaL[1]) && (qt > cutQTL) && (qt < q)) {
        //Treat as Lambda (pos proton, neg pion)
        fillV0Histos<kPr>(posTrack, posTrack.tpcInnerParam(), posTrack.tpcExpSignalDiffPr(), posTrack.tpcNSigmaPr());
        fillV0Histos<kPi>(negTrack, negTrack.tpcInnerParam(), negTrack.tpcExpSignalDiffPi(), negTrack.tpcNSigmaPi());
      }

      // Check for antilambda
      q = cutAPL[0] * sqrt(abs(1 - ((alpha - cutAPL[1]) * (alpha - cutAPL[1])) / (cutAPL[2] * cutAPL[2])));
      if ((alpha > cutAlphaAL[0]) && (alpha < cutAlphaAL[1]) && (qt < q)) {
        //Treat as antilambda (pos pion, neg proton)
        fillV0Histos<kPi>(posTrack, posTrack.tpcInnerParam(), posTrack.tpcExpSignalDiffPi(), posTrack.tpcNSigmaPi());
        fillV0Histos<kPr>(negTrack, negTrack.tpcInnerParam(), negTrack.tpcExpSignalDiffPr(), negTrack.tpcNSigmaPr());
      }

    } //for
  }   //process

}; //struct QaTpcV0

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w;
  w.push_back(adaptAnalysisTask<QaTpcTof>(cfgc));
  if (cfgc.options().get<int>("useV0")) {
    w.push_back(adaptAnalysisTask<QaTpcV0>(cfgc));
  }

  return w;
}