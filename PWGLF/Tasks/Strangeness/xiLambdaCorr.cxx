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

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TDatabasePDG.h"
#include <Math/Vector4D.h>

#include <random>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using Collisions = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>::iterator;
using FullV0s = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras>;
using FullCascades = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs>;
using TracksFull = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

struct xiLambdaCorr {
  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  // binning of (anti)lambda mass QA histograms
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {100, o2::constants::physics::MassLambda - 0.05f, o2::constants::physics::MassLambda + 0.05f}, "binning for the lambda invariant-mass"};
  ConfigurableAxis massXiAxis{"massXiAxis", {100, o2::constants::physics::MassXiMinus - 0.05f, o2::constants::physics::MassXiMinus + 0.05f}, "binning for the Xi invariant-mass"};
  ConfigurableAxis massXiLambdaAxis{"massXiLambdaAxis", {200, o2::constants::physics::MassXiMinus + o2::constants::physics::MassLambda, o2::constants::physics::MassXiMinus + o2::constants::physics::MassLambda + 0.1f}, "binning for the Xi+Lambda invariant-mass"};
  ConfigurableAxis cosPAxis{"cosPAxis", {10, 0.99f, 1.f}, "binning for the cosine of the pointing angle"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.9f, "maximum eta"};
  ConfigurableAxis momAxis{"momAxisFine", {50, 0.f, 10.f}, "momentum axis binning"};
  ConfigurableAxis mixTypeAxis{"mixTypeAxis", {4, -0.5f, 3.5f}, "mixing type axis"}; // xi - lambda , xi - anti-lambda, anti-xi - lambda, anti-xi - anti-lambda

  Configurable<float> cascPtMin{"cascPtMin", 1.f, "minimum (anti)casc pT (GeV/c)"};
  Configurable<float> cascPtMax{"cascPtMax", 4.f, "maximum (anti)casc pT (GeV/c)"};

  Configurable<float> minNTPCClus{"minNTPCClus", 100, "Minimum number of TPC clusters"};
  Configurable<double> minCascCosPA{"minCascCosPA", 0.99f, "Minimum cosine of the pointing angle of the cascade"};
  Configurable<double> minV0CosPA{"minV0CosPA", 0.97f, "Minimum cosine of the pointing angle of the V0"};

  Configurable<float> nSigmaTPCCut{"nSigmaTPCCut", 3.f, "Number of sigmas for the TPC PID"};
  Configurable<float> dcaBachToPV{"dcaBachToPV", 0.05f, "DCA of the bachelor to the primary vertex"};
  Configurable<float> dcaV0Bach{"dcaV0Bach", 1.f, "DCA between the V0 daughters"};

  Configurable<float> lambdaPtMin{"lambdaPtMin", 0.5f, "minimum (anti)lambda pT (GeV/c)"};
  Configurable<float> lambdaPtMax{"lambdaPtMax", 4.f, "maximum (anti)lambda pT (GeV/c)"};
  Configurable<float> dcaLambdaDauToPV{"dcaLambdaDauToPV", 0.05f, "DCA of the lambda daughter to the primary vertex"};
  Configurable<float> minLambdaCosPA{"minLambdaCosPA", 0.99f, "Minimum cosine of the pointing angle of the lambda"};

  Configurable<float> mLambdaWindow{"mLambdaWindow", 0.02f, "mLambdaWindow"};
  Configurable<float> mXiWindow{"mXiWindow", 0.02f, "mXiWindow"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  template <class T>
  bool selectTrack(T const& track)
  {
    if (track.tpcNClsFound() < minNTPCClus) {
      return false;
    }
    return true;
  }

  void init(o2::framework::InitContext&)
  {
    // event QA
    histos.add<TH1>("QA/zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    histos.add<TH2>("QA/massLambda", ";#it{p}_{T} (GeV/#it{c});#it{m}_{#Lambda} (GeV/#it{c}^{2})", HistType::kTH2F, {momAxis, massLambdaAxis});
    histos.add<TH2>("QA/massXi", ";#it{p}_{T} (GeV/#it{c});#it{m}_{#Xi} (GeV/#it{c}^{2})", HistType::kTH2F, {momAxis, massXiAxis});
    histos.add("xiMinusLambda", "", {HistType::kTHnSparseF, {massXiLambdaAxis, momAxis, massXiAxis, massLambdaAxis, cosPAxis, cosPAxis, mixTypeAxis}});
  }

  template <class C, class T>
  bool isSelectedCasc(C const& collision, T const&, FullCascades::iterator const& casc)
  {

    if (std::abs(casc.positiveeta()) > 0.9 || std::abs(casc.negativeeta()) > 0.9 || std::abs(casc.bacheloreta()) > 0.9) {
      return false;
    }

    auto bachelor = casc.bachTrackExtra_as<T>();
    auto posDau = casc.posTrackExtra_as<T>();
    auto negDau = casc.negTrackExtra_as<T>();

    if (!selectTrack(bachelor) || !selectTrack(posDau) || !selectTrack(negDau)) {
      return false;
    }
    if (casc.sign() > 0) {
      if (TMath::Abs(posDau.tpcNSigmaPi()) > nSigmaTPCCut || TMath::Abs(negDau.tpcNSigmaPr()) > nSigmaTPCCut) {
        return false;
      }
    } else if (casc.sign() < 0) {
      if (TMath::Abs(negDau.tpcNSigmaPi()) > nSigmaTPCCut || TMath::Abs(posDau.tpcNSigmaPr()) > nSigmaTPCCut) {
        return false;
      }
    }
    if (TMath::Abs(casc.dcabachtopv()) < dcaBachToPV) {
      return false;
    }
    if (TMath::Abs(casc.dcacascdaughters()) > dcaV0Bach) {
      return false;
    }
    if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < minCascCosPA) {
      return false;
    }
    if (TMath::Abs(casc.eta()) > etaMax) {
      return false;
    }
    // mass cuts
    bool massInWindow = false;
    if (casc.mXi() > o2::constants::physics::MassXiMinus - mXiWindow && casc.mXi() < o2::constants::physics::MassXiMinus + mXiWindow) {
      massInWindow = true;
    }

    if (!massInWindow) {
      return false;
    }
    histos.fill(HIST("QA/massXi"), casc.pt(), casc.mXi());
    return true;
  };

  template <class T>
  bool isSelectedLambda(T const&, FullV0s::iterator const& v0)
  {
    auto posDau = v0.posTrackExtra_as<T>();
    auto negDau = v0.negTrackExtra_as<T>();

    if (std::abs(v0.positiveeta()) > 0.9 || std::abs(v0.negativeeta()) > 0.9) {
      return false;
    }

    if (!selectTrack(posDau) || !selectTrack(negDau)) {
      return false;
    }
    if (v0.alpha() > 0) {
      if (TMath::Abs(posDau.tpcNSigmaPr()) > nSigmaTPCCut || TMath::Abs(negDau.tpcNSigmaPi()) > nSigmaTPCCut) {
        return false;
      }
    } else {
      if (TMath::Abs(posDau.tpcNSigmaPi()) > nSigmaTPCCut || TMath::Abs(negDau.tpcNSigmaPr()) > nSigmaTPCCut) {
        return false;
      }
    }

    if (TMath::Abs(v0.dcapostopv()) < dcaLambdaDauToPV || TMath::Abs(v0.dcanegtopv()) < dcaLambdaDauToPV) {
      return false;
    }

    if (v0.v0cosPA() < minLambdaCosPA) {
      return false;
    }
    if (TMath::Abs(v0.eta()) > etaMax) {
      return false;
    }
    // mass cuts
    bool massInWindow = false;
    float massLambda = v0.alpha() > 0 ? v0.mLambda() : v0.mAntiLambda();
    if (v0.mLambda() > o2::constants::physics::MassLambda - mLambdaWindow && v0.mLambda() < o2::constants::physics::MassLambda + mLambdaWindow) {
      massInWindow = true;
    }
    if (!massInWindow) {
      return false;
    }
    histos.fill(HIST("QA/massLambda"), v0.pt(), massLambda);
    return true;
  };

  template <class C, class T>
  void fillXiLambda(C const& collision, T const& tracks, FullV0s const& v0s, FullCascades const& cascades)
  {
    for (auto& casc : cascades) {
      if (!isSelectedCasc(collision, tracks, casc)) {
        continue;
      }
      for (auto& v0 : v0s) {
        if (!isSelectedLambda(tracks, v0)) {
          continue;
        }
        if (casc.posTrackExtraId() == v0.posTrackExtraId() || casc.posTrackExtraId() == v0.negTrackExtraId()) {
          continue;
        }

        int mixType = -1;
        if (casc.sign() > 0 && v0.alpha() > 0) {
          mixType = 0; // xi - lambda
        } else if (casc.sign() > 0 && v0.alpha() < 0) {
          mixType = 1; // xi - anti-lambda
        } else if (casc.sign() < 0 && v0.alpha() > 0) {
          mixType = 2; // anti-xi - lambda
        } else if (casc.sign() < 0 && v0.alpha() < 0) {
          mixType = 3; // anti-xi - anti-lambda
        }

        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>> cascMom4D(casc.px(), casc.py(), casc.pz(), o2::constants::physics::MassXiMinus);
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>> lambdaMom4D{v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassLambda};
        auto xiLambdaMom4D = cascMom4D + lambdaMom4D;
        float massLambda = v0.alpha() > 0 ? v0.mLambda() : v0.mAntiLambda();
        float massXi = casc.mXi();
        float cosPAxi = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
        float cosPAlambda = v0.v0cosPA();
        histos.fill(HIST("xiMinusLambda"), xiLambdaMom4D.M(), xiLambdaMom4D.Pt(), massXi, massLambda, cosPAxi, cosPAlambda, mixType);
      }
    }
  };

  void processData(Collisions const& collision, TracksFull const& tracks, FullV0s const& v0s, FullCascades const& cascades)
  {

    if (!collision.sel8())
      return;

    if (std::abs(collision.posZ()) > zVtxMax)
      return;

    if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return;

    histos.fill(HIST("QA/zVtx"), collision.posZ());
    fillXiLambda(collision, tracks, v0s, cascades);
  }
  PROCESS_SWITCH(xiLambdaCorr, processData, "Process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<xiLambdaCorr>(cfgc)};
}
