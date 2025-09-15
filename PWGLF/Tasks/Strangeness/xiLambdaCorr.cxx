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

#include "PWGLF/DataModel/LFDoubleCascTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TDatabasePDG.h"
#include <Math/Vector4D.h>

#include <random>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults>::iterator;
using FullCascades = aod::CascDataExt;
using FullV0s = aod::V0Datas;
using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;
struct xiLambdaCorr {
  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};
  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  // binning of (anti)lambda mass QA histograms
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {200, o2::constants::physics::MassLambda - 0.05f, o2::constants::physics::MassLambda + 0.05f}, "binning for the lambda invariant-mass"};
  ConfigurableAxis massXiAxis{"massXiAxis", {200, o2::constants::physics::MassXiMinus - 0.05f, o2::constants::physics::MassXiMinus + 0.05f}, "binning for the Xi invariant-mass"};
  ConfigurableAxis massXiLambdaAxis{"massXiLambdaAxis", {200, o2::constants::physics::MassXiMinus + o2::constants::physics::MassLambda, o2::constants::physics::MassXiMinus + o2::constants::physics::MassLambda + 0.1f}, "binning for the Xi+Lambda invariant-mass"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.9f, "maximum eta"};
  ConfigurableAxis momAxis{"momAxisFine", {5.e2, 0.f, 5.f}, "momentum axis binning"};
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
    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
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
    histos.add("xiMinusLambda", "", {HistType::kTHnSparseF, {massXiLambdaAxis, momAxis, mixTypeAxis}});
  }

  template <class C, class T>
  bool isSelectedCasc(C const& collision, T const&, FullCascades::iterator const& casc)
  {

    auto bachelor = casc.bachelor_as<T>();
    auto posDau = casc.posTrack_as<T>();
    auto negDau = casc.negTrack_as<T>();

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
    auto posDau = v0.posTrack_as<T>();
    auto negDau = v0.negTrack_as<T>();

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
        if (casc.posTrackId() == v0.posTrackId() || casc.posTrackId() == v0.negTrackId()) {
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
        histos.fill(HIST("xiMinusLambda"), xiLambdaMom4D.M(), xiLambdaMom4D.Pt(), mixType);
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<xiLambdaCorr>(cfgc)};
}
