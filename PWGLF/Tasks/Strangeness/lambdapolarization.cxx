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

/// \author Junlee Kim (jikim1290@gmail.com)

#include <cmath>
#include <array>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TVector2.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include <TMath.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"

#include "CommonConstants/PhysicsConstants.h"

#include "ReconstructionDataFormats/Track.h"

#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "CCDB/BasicCCDBManager.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct lambdapolarization {

  //  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults, aod::Qvectors>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using V0TrackCandidate = aod::V0Datas;

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgv0radiusMin{"cfgv0radiusMin", 1.2, "minimum decay radius"};
  Configurable<float> cfgDCAPosToPVMin{"cfgDCAPosToPVMin", 0.05, "minimum DCA to PV for positive track"};
  Configurable<float> cfgDCANegToPVMin{"cfgDCANegToPVMin", 0.2, "minimum DCA to PV for negative track"};
  Configurable<float> cfgv0CosPA{"cfgv0CosPA", 0.995, "minimum v0 cosine"};
  Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};

  Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
  Configurable<float> cfgV0EtaMin{"cfgV0EtaMin", -0.5, "maximum rapidity"};
  Configurable<float> cfgV0EtaMax{"cfgV0EtaMax", 0.5, "maximum rapidity"};
  Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};

  Configurable<bool> cfgQAv0{"cfgQAv0", false, "QA plot"};

  Configurable<int> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 70, "minimum fired crossed rows"};
  Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 5, "proton nsigma for TPC"};
  Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 5, "pion nsigma for TPC"};
  Configurable<float> cfgDaughEtaMin{"cfgDaughEtaMin", -0.8, "minimum daughter eta"};
  Configurable<float> cfgDaughEtaMax{"cfgDaughEtaMax", 0.8, "maximum daughter eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.6, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};

  Configurable<int> cfgnMods{"cfgnMods", 1, "The number of modulations of interest starting from 2"};

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;

  double angle;
  double relphi;

  void init(o2::framework::InitContext&)
  {
    if (cfgnMods > 1)
      LOGF(fatal, "multiple harmonics not implemented yet"); // FIXME: will be updated after Qvector task updates

    AxisSpec massAxis = {100, 1.0, 1.2};
    AxisSpec ptAxis = {100, 0.0, 10.0};
    AxisSpec cosAxis = {100, -1.0, 1.0};
    AxisSpec centAxis = {80, 0.0, 80.0};
    AxisSpec epAxis = {6, 0.0, 2.0 * constants::math::PI};

    AxisSpec pidAxis = {100, -10, 10};
    for (auto i = 2; i < cfgnMods + 2; i++) {
      histos.add(Form("h_lambda_cos_psi%d", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("h_alambda_cos_psi%d", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("h_lambda_cos2_psi%d", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("h_alambda_cos2_psi%d", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
    }

    if (cfgQAv0) {
      histos.add("QA/nsigma_tpc_pt_ppr", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_ppi", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_mpr", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_mpi", "", {HistType::kTH2F, {ptAxis, pidAxis}});
    }

    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
  }

  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  ROOT::Math::PxPyPzMVector ProtonVec, PionVec, LambdaVec, ProtonBoostedVec, PionBoostedVec;

  template <typename TCollision>
  bool eventSelected(TCollision collision)
  {
    if (!collision.sel8()) {
      return 0;
    }
    auto centrality = collision.centFT0C();
    auto multNTracksPV = collision.multNTracksPV();

    if (multNTracksPV < fMultPVCutLow->Eval(centrality)) {
      return 0;
    }

    if (multNTracksPV > fMultPVCutHigh->Eval(centrality)) {
      return 0;
    }

    if ((!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return 0;
    }

    return 1;
  } // event selection

  template <typename TCollision, typename V0>
  bool SelectionV0(TCollision const& collision, V0 const& candidate)
  {
    if (candidate.v0radius() < cfgv0radiusMin)
      return false;
    if (candidate.dcapostopv() < cfgDCAPosToPVMin)
      return false;
    if (candidate.dcanegtopv() < cfgDCANegToPVMin)
      return false;
    if (candidate.v0cosPA() < cfgv0CosPA)
      return false;
    if (candidate.dcaV0daughters() > cfgDCAV0Dau)
      return false;
    if (candidate.pt() < cfgV0PtMin)
      return false;
    if (candidate.eta() < cfgV0EtaMin)
      return false;
    if (candidate.eta() > cfgV0EtaMax)
      return false;
    if (candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda > cfgV0LifeTime)
      return false;

    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, int pid) // pid 0: proton, pid 1: pion
  {
    if (track.tpcNClsFound() < cfgDaughTPCnclsMin)
      return false;
    if (pid == 0 && std::abs(track.tpcNSigmaPr()) > cfgDaughPIDCutsTPCPr)
      return false;
    if (pid == 1 && std::abs(track.tpcNSigmaPi()) > cfgDaughPIDCutsTPCPi)
      return false;
    if (track.eta() > cfgDaughEtaMax)
      return false;
    if (track.eta() < cfgDaughEtaMin)
      return false;
    if (pid == 0 && track.pt() < cfgDaughPrPt)
      return false;
    if (pid == 1 && track.pt() < cfgDaughPiPt)
      return false;

    return true;
  }

  template <typename TCollision, typename V0>
  void FillHistograms(TCollision const& collision, V0 const& V0s, int nmode)
  {
    for (auto& v0 : V0s) {
      auto postrack = v0.template posTrack_as<TrackCandidates>();
      auto negtrack = v0.template negTrack_as<TrackCandidates>();

      double nTPCSigmaPosPr = postrack.tpcNSigmaPr();
      double nTPCSigmaNegPi = negtrack.tpcNSigmaPi();

      double nTPCSigmaNegPr = negtrack.tpcNSigmaPr();
      double nTPCSigmaPosPi = postrack.tpcNSigmaPi();

      if (cfgQAv0) {
        histos.fill(HIST("QA/nsigma_tpc_pt_ppr"), postrack.pt(), nTPCSigmaPosPr);
        histos.fill(HIST("QA/nsigma_tpc_pt_ppi"), postrack.pt(), nTPCSigmaPosPi);

        histos.fill(HIST("QA/nsigma_tpc_pt_mpr"), negtrack.pt(), nTPCSigmaNegPr);
        histos.fill(HIST("QA/nsigma_tpc_pt_mpi"), negtrack.pt(), nTPCSigmaNegPi);
      }

      int LambdaTag = 0;
      int aLambdaTag = 0;

      if (isSelectedV0Daughter(postrack, 0) && isSelectedV0Daughter(negtrack, 1)) {
        LambdaTag = 1;
      }
      if (isSelectedV0Daughter(negtrack, 0) && isSelectedV0Daughter(postrack, 1)) {
        aLambdaTag = 1;
      }

      if (!LambdaTag && !aLambdaTag)
        continue;

      if (!SelectionV0(collision, v0))
        continue;

      if (LambdaTag) {
        ProtonVec = ROOT::Math::PxPyPzMVector(postrack.px(), postrack.py(), postrack.pz(), massPr);
        PionVec = ROOT::Math::PxPyPzMVector(negtrack.px(), negtrack.py(), negtrack.pz(), massPi);
      }
      if (aLambdaTag) {
        ProtonVec = ROOT::Math::PxPyPzMVector(negtrack.px(), negtrack.py(), negtrack.pz(), massPr);
        PionVec = ROOT::Math::PxPyPzMVector(postrack.px(), postrack.py(), postrack.pz(), massPi);
      }
      LambdaVec = ProtonVec + PionVec;

      ROOT::Math::Boost boost{LambdaVec.BoostToCM()};
      ProtonBoostedVec = boost(ProtonVec);

      angle = ProtonBoostedVec.Pz() / ProtonBoostedVec.P();
      relphi = TVector2::Phi_0_2pi(static_cast<float>(nmode) * (LambdaVec.Phi() - TMath::ATan2(collision.qvecIm()[3 + (nmode - 2) * 24], collision.qvecRe()[3 + (nmode - 2) * 24])));

      if (nmode == 2) { ////////////
        if (LambdaTag) {
          histos.fill(HIST("h_lambda_cos_psi2"), v0.mLambda(), v0.pt(), angle, collision.centFT0C(), relphi);
          histos.fill(HIST("h_lambda_cos2_psi2"), v0.mLambda(), v0.pt(), angle * angle, collision.centFT0C(), relphi);
        }
        if (aLambdaTag) {
          histos.fill(HIST("h_alambda_cos_psi2"), v0.mLambda(), v0.pt(), angle, collision.centFT0C(), relphi);
          histos.fill(HIST("h_alambda_cos2_psi2"), v0.mLambda(), v0.pt(), angle * angle, collision.centFT0C(), relphi);
        }
      } else if (nmode == 3) {
        if (LambdaTag) {
          histos.fill(HIST("h_lambda_cos_psi3"), v0.mLambda(), v0.pt(), angle, collision.centFT0C(), relphi);
          histos.fill(HIST("h_lambda_cos2_psi3"), v0.mLambda(), v0.pt(), angle * angle, collision.centFT0C(), relphi);
        }
        if (aLambdaTag) {
          histos.fill(HIST("h_alambda_cos_psi3"), v0.mLambda(), v0.pt(), angle, collision.centFT0C(), relphi);
          histos.fill(HIST("h_alambda_cos2_psi3"), v0.mLambda(), v0.pt(), angle * angle, collision.centFT0C(), relphi);
        }
      } ////////// FIXME: not possible to get histograms using nmode
    }
  }

  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks, aod::V0Datas const& V0s,
                   aod::BCs const&)
  {
    if (!eventSelected(collision)) {
      return;
    }
    for (auto i = 2; i < cfgnMods + 2; i++) {
      FillHistograms(collision, V0s, i);
    }
  }
  PROCESS_SWITCH(lambdapolarization, processData, "Process Event for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdapolarization>(cfgc, TaskName{"lf-lambdapolarization"})};
}
