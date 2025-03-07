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

/// \author Youpeng Su (yousu@cern.ch)
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include <TLorentzVector.h>
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include <TTree.h>
#include <TFile.h>
#include <TMatrixD.h>
#include "PWGLF/DataModel/lambdaJetpolarization.h"
using std::cout;
using std::endl;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct LfMyV0s {
  HistogramRegistry registry{"registry"};
  void init(InitContext const&)
  {
    const AxisSpec axisPx{100, -10, 10, "#px (GeV/c)"};
    const AxisSpec axisPy{100, -10, 10, "#py (GeV/c)"};
    const AxisSpec axisPz{100, -10, 10, "#pz (GeV/c)"};
    const AxisSpec axisPT{200, 0, 50, "#p_{T} (GeV/c)"};
    const AxisSpec axisPhi{100, -3.14, 3.14, "#Phi"};
    const AxisSpec axisMass{100, 0, 2, "Mass(GeV/c^{2})"};

    registry.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 0.9f, 1.2f}}});
    registry.add("V0pTInLab", "V0pTInLab", kTH1F, {axisPT});
    registry.add("hMassVsPtLambda", "hMassVsPtLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});

    registry.add("V0pxInLab", "V0pxInLab", kTH1F, {axisPx});
    registry.add("V0pyInLab", "V0pyInLab", kTH1F, {axisPy});
    registry.add("V0pzInLab", "V0pzInLab", kTH1F, {axisPz});

    registry.add("V0pxInRest_frame", "V0pxInRest_frame", kTH1F, {axisPx});
    registry.add("V0pyInRest_frame", "V0pyInRest_frame", kTH1F, {axisPy});
    registry.add("V0pzInRest_frame", "V0pzInRest_frame", kTH1F, {axisPz});

    registry.add("JetpxInLab", "JetpxInLab", kTH1F, {axisPx});
    registry.add("JetpyInLab", "JetpyInLab", kTH1F, {axisPy});
    registry.add("JetpzInLab", "JetpzInLab", kTH1F, {axisPz});
    registry.add("JetpTInLab", "JetpTInLab", kTH1F, {axisPT});

    registry.add("LeadingJetpx", "LeadingJetpx", kTH1F, {axisPx});
    registry.add("LeadingJetpy", "LeadingJetpy", kTH1F, {axisPy});
    registry.add("LeadingJetpz", "LeadingJetpz", kTH1F, {axisPz});
    registry.add("LeadingJetpT", "LeadingJetpT", kTH1F, {axisPT});

    registry.add("V0protonpxInLab", "V0protonpxInLab", kTH1F, {axisPx});
    registry.add("V0protonpyInLab", "V0protonpyInLab", kTH1F, {axisPy});
    registry.add("V0protonpzInLab", "V0protonpzInLab", kTH1F, {axisPz});
    registry.add("V0protonphiInLab", "V0protonphiInLab", kTH1F, {axisPhi});

    registry.add("V0protonpxInRest_frame", "V0protonpxInRest_frame", kTH1F, {axisPx});
    registry.add("V0protonpyInRest_frame", "V0protonpyInRest_frame", kTH1F, {axisPy});
    registry.add("V0protonpzInRest_frame", "V0protonpzInRest_frame", kTH1F, {axisPz});
    registry.add("V0protonMassInRest_frame", "V0protonMassInRest_frame", kTH1F, {axisMass});
    registry.add("V0protonphiInRest_frame", "V0protonphiInRest_frame", kTH1F, {axisPhi});

    registry.add("V0protonpxInJetV0frame", "V0protonpxInJetV0frame", kTH1F, {axisPx});
    registry.add("V0protonpyInJetV0frame", "V0protonpyInJetV0frame", kTH1F, {axisPy});
    registry.add("V0protonpzInJetV0frame", "V0protonpzInJetV0frame", kTH1F, {axisPz});
    registry.add("V0protonphiInJetV0frame", "V0protonphiInJetV0frame", kTH1F, {axisPhi});

    registry.add("V0LambdapxInJetV0frame", "V0LambdapxInJetV0frame", kTH1F, {axisPx});
    registry.add("V0LambdapyInJetV0frame", "V0LambdapyInJetV0frame", kTH1F, {axisPy});
    registry.add("V0LambdapzInJetV0frame", "V0LambdapzInJetV0frame", kTH1F, {axisPz});

    registry.add("hLambdamassandSinPhi", "hLambdamassandSinPhi", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registry.add("profile", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
  }
  double massPr = o2::constants::physics::MassProton;
  double massLambda = o2::constants::physics::MassLambda;

  TMatrixD LorentzTransInV0frame(double ELambda, double Lambdapx, double Lambdapy, double Lambdapz)
  {
    double PLambda = sqrt(Lambdapx * Lambdapx + Lambdapy * Lambdapy + Lambdapz * Lambdapz);
    double LambdaMass = sqrt(ELambda * ELambda - PLambda * PLambda);
    double Alpha = 1 / (LambdaMass * (ELambda + LambdaMass));
    TMatrixD matrixLabToLambda(4, 4);
    matrixLabToLambda(0, 0) = ELambda / LambdaMass;
    matrixLabToLambda(0, 1) = -Lambdapx / LambdaMass;
    matrixLabToLambda(0, 2) = -Lambdapy / LambdaMass;
    matrixLabToLambda(0, 3) = -Lambdapz / LambdaMass;
    matrixLabToLambda(1, 0) = -Lambdapx / LambdaMass;
    matrixLabToLambda(1, 1) = 1 + Alpha * Lambdapx * Lambdapx;
    matrixLabToLambda(1, 2) = Alpha * Lambdapx * Lambdapy;
    matrixLabToLambda(1, 3) = Alpha * Lambdapx * Lambdapz;
    matrixLabToLambda(2, 0) = -Lambdapy / LambdaMass;
    matrixLabToLambda(2, 1) = Alpha * Lambdapy * Lambdapx;
    matrixLabToLambda(2, 2) = 1 + Alpha * Lambdapy * Lambdapy;
    matrixLabToLambda(2, 3) = Alpha * Lambdapy * Lambdapz;
    matrixLabToLambda(3, 0) = -Lambdapz / LambdaMass;
    matrixLabToLambda(3, 1) = Alpha * Lambdapz * Lambdapx;
    matrixLabToLambda(3, 2) = Alpha * Lambdapz * Lambdapy;
    matrixLabToLambda(3, 3) = 1 + Alpha * Lambdapz * Lambdapz;
    return matrixLabToLambda;
  }
  TMatrixD MyTMatrixTranslationToJet(double Jetpx, double Jetpy, double Jetpz, double Lambdapx, double Lambdapy, double Lambdapz)
  {
    TVector3 UnitX(1.0, 0.0, 0.0);
    TVector3 UnitY(0.0, 1.0, 0.0);
    TVector3 UnitZ(0.0, 0.0, 1.0);
    TVector3 JetP(Jetpx, Jetpy, Jetpz);
    TVector3 V0LambdaP(Lambdapx, Lambdapy, Lambdapz);
    TVector3 vortex_y = (JetP.Cross(V0LambdaP));

    TVector3 z_hat = JetP.Unit();
    TVector3 y_hat = vortex_y.Unit();
    TVector3 x_hat1 = y_hat.Cross(z_hat);
    TVector3 x_hat = x_hat1.Unit();

    TMatrixD matrixLabToJet(4, 4);
    matrixLabToJet(0, 0) = 1;
    matrixLabToJet(0, 1) = 0.0;
    matrixLabToJet(0, 2) = 0.0;
    matrixLabToJet(0, 3) = 0.0;
    matrixLabToJet(1, 0) = 0.0;
    matrixLabToJet(1, 1) = x_hat.X();
    matrixLabToJet(1, 2) = x_hat.Y();
    matrixLabToJet(1, 3) = x_hat.Z();
    matrixLabToJet(2, 0) = 0.0;
    matrixLabToJet(2, 1) = y_hat.X();
    matrixLabToJet(2, 2) = y_hat.Y();
    matrixLabToJet(2, 3) = y_hat.Z();
    matrixLabToJet(3, 0) = 0.0;
    matrixLabToJet(3, 1) = z_hat.X();
    matrixLabToJet(3, 2) = z_hat.Y();
    matrixLabToJet(3, 3) = z_hat.Z();
    return matrixLabToJet;
  }
  // aod::MyCollision const& collision
  void processJetV0Analysis(aod::MyTable const& myv0s, aod::MyTableJet const& myJets)
  {
    for (auto& candidate : myv0s) {
      registry.fill(HIST("hMassLambda"), candidate.v0Lambdamass());
      registry.fill(HIST("V0pTInLab"), candidate.v0pt());
      registry.fill(HIST("hMassVsPtLambda"), candidate.v0pt(), candidate.v0Lambdamass());
      registry.fill(HIST("V0pxInLab"), candidate.v0px());
      registry.fill(HIST("V0pyInLab"), candidate.v0py());
      registry.fill(HIST("V0pzInLab"), candidate.v0pz());
      registry.fill(HIST("V0protonpxInLab"), candidate.v0protonpx());
      registry.fill(HIST("V0protonpyInLab"), candidate.v0protonpy());
      registry.fill(HIST("V0protonpzInLab"), candidate.v0protonpz());
      double protonsinPhiInLab = candidate.v0protonpy() / sqrt(candidate.v0protonpx() * candidate.v0protonpx() + candidate.v0protonpy() * candidate.v0protonpy());
      registry.fill(HIST("V0protonphiInLab"), protonsinPhiInLab);
      double PLambda = sqrt(candidate.v0px() * candidate.v0px() + candidate.v0py() * candidate.v0py() + candidate.v0pz() * candidate.v0pz());
      double ELambda = sqrt(candidate.v0Lambdamass() * candidate.v0Lambdamass() + PLambda * PLambda);
      TMatrixD pLabV0(4, 1);
      pLabV0(0, 0) = ELambda;
      pLabV0(1, 0) = candidate.v0px();
      pLabV0(2, 0) = candidate.v0py();
      pLabV0(3, 0) = candidate.v0pz();
      TMatrixD V0InV0(4, 1);
      V0InV0 = LorentzTransInV0frame(ELambda, candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabV0;
      registry.fill(HIST("V0pxInRest_frame"), V0InV0(1, 0));
      registry.fill(HIST("V0pyInRest_frame"), V0InV0(2, 0));
      registry.fill(HIST("V0pzInRest_frame"), V0InV0(3, 0));
    }
    for (auto& candidate : myv0s) {
      double PLambda = sqrt(candidate.v0px() * candidate.v0px() + candidate.v0py() * candidate.v0py() + candidate.v0pz() * candidate.v0pz());
      double ELambda = sqrt(candidate.v0Lambdamass() * candidate.v0Lambdamass() + PLambda * PLambda);
      TMatrixD pLabproton(4, 1);
      double protonE = sqrt(massPr * massPr + candidate.v0protonpx() * candidate.v0protonpx() + candidate.v0protonpy() * candidate.v0protonpy() + candidate.v0protonpz() * candidate.v0protonpz());
      pLabproton(0, 0) = protonE;
      pLabproton(1, 0) = candidate.v0protonpx();
      pLabproton(2, 0) = candidate.v0protonpy();
      pLabproton(3, 0) = candidate.v0protonpz();
      TMatrixD protonInV0(4, 1);
      protonInV0 = LorentzTransInV0frame(ELambda, candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabproton;
      double protonMassInV0 = sqrt(protonInV0(0, 0) * protonInV0(0, 0) - protonInV0(1, 0) * protonInV0(1, 0) - protonInV0(2, 0) * protonInV0(2, 0) - protonInV0(3, 0) * protonInV0(3, 0));
      registry.fill(HIST("V0protonMassInRest_frame"), protonMassInV0);
      registry.fill(HIST("V0protonpxInRest_frame"), protonInV0(1, 0));
      registry.fill(HIST("V0protonpyInRest_frame"), protonInV0(2, 0));
      registry.fill(HIST("V0protonpzInRest_frame"), protonInV0(3, 0));
      double protonsinPhiInV0frame = protonInV0(2, 0) / sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0));
      registry.fill(HIST("V0protonphiInRest_frame"), protonsinPhiInV0frame);
    }

    for (auto& Jet : myJets) {
      registry.fill(HIST("JetpxInLab"), Jet.jetpx());
      registry.fill(HIST("JetpyInLab"), Jet.jetpy());
      registry.fill(HIST("JetpzInLab"), Jet.jetpz());
      registry.fill(HIST("JetpTInLab"), Jet.jetpt());
    }
  }
  PROCESS_SWITCH(LfMyV0s, processJetV0Analysis, "processJetV0Analysis", true);
  void processLeadingJetV0Analysis(aod::MyTable const& myv0s, aod::MyTableLeadingJet const& myleadingJets)
  {
    for (auto& LeadingJet : myleadingJets) {
      int V0Numbers = 0;
      double protonsinPhiInJetV0frame = 0;
      for (auto& candidate : myv0s) {
        if (candidate.mycollisionv0() == LeadingJet.mycollisionleadingjet()) {
          V0Numbers = V0Numbers + 1;
          double PLambda = sqrt(candidate.v0px() * candidate.v0px() + candidate.v0py() * candidate.v0py() + candidate.v0pz() * candidate.v0pz());
          double ELambda = sqrt(candidate.v0Lambdamass() * candidate.v0Lambdamass() + PLambda * PLambda);
          double protonE = sqrt(massPr * massPr + candidate.v0protonpx() * candidate.v0protonpx() + candidate.v0protonpy() * candidate.v0protonpy() + candidate.v0protonpz() * candidate.v0protonpz());

          TMatrixD pLabV0(4, 1);
          pLabV0(0, 0) = ELambda;
          pLabV0(1, 0) = candidate.v0px();
          pLabV0(2, 0) = candidate.v0py();
          pLabV0(3, 0) = candidate.v0pz();

          TMatrixD lambdaInJet(4, 1);
          lambdaInJet = MyTMatrixTranslationToJet(LeadingJet.leadingjetpx(), LeadingJet.leadingjetpy(), LeadingJet.leadingjetpz(), candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabV0;

          TMatrixD lambdaInJetV0(4, 1);
          lambdaInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(LeadingJet.leadingjetpx(), LeadingJet.leadingjetpy(), LeadingJet.leadingjetpz(), candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabV0;
          registry.fill(HIST("V0LambdapxInJetV0frame"), lambdaInJetV0(1, 0));
          registry.fill(HIST("V0LambdapyInJetV0frame"), lambdaInJetV0(2, 0));
          registry.fill(HIST("V0LambdapzInJetV0frame"), lambdaInJetV0(3, 0));

          TMatrixD pLabproton(4, 1);
          pLabproton(0, 0) = protonE;
          pLabproton(1, 0) = candidate.v0protonpx();
          pLabproton(2, 0) = candidate.v0protonpy();
          pLabproton(3, 0) = candidate.v0protonpz();
          TMatrixD protonInJetV0(4, 1);
          protonInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(LeadingJet.leadingjetpx(), LeadingJet.leadingjetpy(), LeadingJet.leadingjetpz(), candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabproton;
          registry.fill(HIST("V0protonpxInJetV0frame"), protonInJetV0(1, 0));
          registry.fill(HIST("V0protonpyInJetV0frame"), protonInJetV0(2, 0));
          registry.fill(HIST("V0protonpzInJetV0frame"), protonInJetV0(3, 0));
          protonsinPhiInJetV0frame = protonsinPhiInJetV0frame + protonInJetV0(2, 0) / sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));
        }
      }
      for (auto& candidate : myv0s) {
        if (candidate.mycollisionv0() == LeadingJet.mycollisionleadingjet()) {
          registry.fill(HIST("V0protonphiInJetV0frame"), protonsinPhiInJetV0frame / V0Numbers);
          registry.fill(HIST("hLambdamassandSinPhi"), candidate.v0Lambdamass(), protonsinPhiInJetV0frame / V0Numbers);
          registry.fill(HIST("profile"), candidate.v0Lambdamass(), protonsinPhiInJetV0frame / V0Numbers);
        }
      }
    }
    for (auto& LeadingJet : myleadingJets) {
      registry.fill(HIST("LeadingJetpx"), LeadingJet.leadingjetpx());
      registry.fill(HIST("LeadingJetpy"), LeadingJet.leadingjetpy());
      registry.fill(HIST("LeadingJetpz"), LeadingJet.leadingjetpz());
      registry.fill(HIST("LeadingJetpT"), LeadingJet.leadingjetpt());
    }
  }
  PROCESS_SWITCH(LfMyV0s, processLeadingJetV0Analysis, "processLeadingJetV0Analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LfMyV0s>(cfgc),
  };
}
