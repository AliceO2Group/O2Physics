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
/// \file   qaLamMomResolution.cxx
/// \author Carolina Reetz c.reetz@cern.ch
/// \brief  QA task to study momentum resolution of Lambda daughter tracks

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DPG/Tasks/AOTTrack/V0Cascades/qaLamMomResolution.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using V0DatasLabeled = soa::Join<aod::V0Datas, aod::V0Covs, aod::V0DauCovs, aod::McV0Labels>;
using MyTracks = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::McTrackLabels>;

struct qaLamMomResolution {

  Produces<aod::LamDaughters> lamdaughters;

  Configurable<bool> collSelection{"collSelection", true, "collision selection"};

  HistogramRegistry hist{"Histograms"};

  int LambdaPDG = 3122;
  int AntiLambdaPDG = -3122;
  int ProtonPDG = 2212;
  int AntiProtonPDG = -2212;
  int PosPionPDG = 211;
  int NegPionPDG = -211;

  float massLambda = -1.0f;
  float radiusLambda = -1.0f;
  float ptLambda = -1.0f;
  float etaProton = -1.0f, etaPion = -1.0f;
  int tpcNClsProton = 0, tpcNClsPion = 0;
  int chargeProton = 0, chargePion = 0;
  // daughter momenta
  std::array<float, 3> momProtonRecIU;
  std::array<float, 3> momPionRecIU;
  std::array<float, 3> momProtonRecIUErr;
  std::array<float, 3> momPionRecIUErr;
  std::array<float, 3> momProtonRec;
  std::array<float, 3> momPionRec;
  std::array<float, 3> momProtonRecErr;
  std::array<float, 3> momPionRecErr;
  float sigma1PtProtonIU = -1.0f, sigma1PtPionIU = -1.0f;
  // daughter DCA
  std::array<float, 2> DCAProtonRec;    // 0: xy, 1: z
  std::array<float, 2> DCAPionRec;      // 0: xy, 1: z
  std::array<float, 2> DCAProtonRecErr; // 0: xy, 1: z
  std::array<float, 2> DCAPionRecErr;   // 0: xy, 1: z
  // MC info
  std::array<float, 3> momProtonGen;
  std::array<float, 3> momPionGen;
  std::array<float, 2> DCAProtonGen; // 0: xy, 1: z
  std::array<float, 2> DCAPionGen;   // 0: xy, 1: z

  void init(InitContext const&)
  {
    auto hEventSelectionFlow = hist.add<TH1>("hEventSelectionFlow", "Event selection flow", kTH1F, {{2, 0.5f, 2.5f}});
    hEventSelectionFlow->GetXaxis()->SetBinLabel(hEventSelectionFlow->FindBin(1), "Sel8");
    hEventSelectionFlow->GetXaxis()->SetBinLabel(hEventSelectionFlow->FindBin(2), "|Vtx_{z}|<10cm");
  }

  template <typename TCollision>
  void fillTable(TCollision const& collision)
  {
    lamdaughters(collision.globalIndex(),
                 massLambda, radiusLambda, ptLambda,
                 chargeProton, chargePion,
                 etaProton, etaPion,
                 tpcNClsProton, tpcNClsPion,
                 momProtonRec[0], momProtonRec[1], momProtonRec[2],
                 momProtonRecErr[0], momProtonRecErr[1], momProtonRecErr[2],
                 momPionRec[0], momPionRec[1], momPionRec[2],
                 momPionRecErr[0], momPionRecErr[1], momPionRecErr[2],
                 momProtonRecIU[0], momProtonRecIU[1], momProtonRecIU[2],
                 momProtonRecIUErr[0], momProtonRecIUErr[1], momProtonRecIUErr[2],
                 momPionRecIU[0], momPionRecIU[1], momPionRecIU[2],
                 momPionRecIUErr[0], momPionRecIUErr[1], momPionRecIUErr[2],
                 momProtonGen[0], momProtonGen[1], momProtonGen[2],
                 momPionGen[0], momPionGen[1], momPionGen[2],
                 sigma1PtProtonIU, sigma1PtPionIU,
                 DCAProtonRec[0], DCAProtonRec[1],
                 DCAProtonRecErr[0], DCAProtonRecErr[1],
                 DCAPionRec[0], DCAPionRec[1],
                 DCAPionRecErr[0], DCAPionRecErr[1]);
  }

  template <typename TProton, typename TPion>
  void getTrackInfo(TProton const& protonTrackIU, TPion const& pionTrackIU)
  {
    // daughter momenta at IU
    o2::track::TrackParCov protonTrackParCov, pionTrackParCov;
    std::array<float, 21> protoncv, pioncv;
    protonTrackParCov = getTrackParCov(protonTrackIU);
    pionTrackParCov = getTrackParCov(pionTrackIU);
    protonTrackParCov.getCovXYZPxPyPzGlo(protoncv);
    pionTrackParCov.getCovXYZPxPyPzGlo(pioncv);
    // proton
    momProtonRecIU[0] = protonTrackIU.px();
    momProtonRecIU[1] = protonTrackIU.py();
    momProtonRecIU[2] = protonTrackIU.pz();
    momProtonRecIUErr[0] = sqrt(protoncv[9]);
    momProtonRecIUErr[1] = sqrt(protoncv[14]);
    momProtonRecIUErr[2] = sqrt(protoncv[20]);
    sigma1PtProtonIU = protonTrackIU.sigma1Pt();
    // pion
    momPionRecIU[0] = pionTrackIU.px();
    momPionRecIU[1] = pionTrackIU.py();
    momPionRecIU[2] = pionTrackIU.pz();
    momPionRecIUErr[0] = sqrt(pioncv[9]);
    momPionRecIUErr[1] = sqrt(pioncv[14]);
    momPionRecIUErr[2] = sqrt(pioncv[20]);
    sigma1PtPionIU = pionTrackIU.sigma1Pt();

    // daughter DCA
    DCAProtonRec[0] = protonTrackIU.dcaXY();
    DCAProtonRec[1] = protonTrackIU.dcaZ();
    DCAProtonRecErr[0] = sqrt(protonTrackIU.sigmaDcaXY2());
    DCAProtonRecErr[1] = sqrt(protonTrackIU.sigmaDcaZ2());
    DCAPionRec[0] = pionTrackIU.dcaXY();
    DCAPionRec[1] = pionTrackIU.dcaZ();
    DCAPionRecErr[0] = sqrt(pionTrackIU.sigmaDcaXY2());
    DCAPionRecErr[1] = sqrt(pionTrackIU.sigmaDcaZ2());

    // daughter charges, eta, nTPCclusters
    chargeProton = protonTrackIU.sign();
    chargePion = pionTrackIU.sign();
    etaProton = protonTrackIU.eta();
    etaPion = pionTrackIU.eta();
    tpcNClsProton = protonTrackIU.tpcNClsFound();
    tpcNClsPion = pionTrackIU.tpcNClsFound();
  }

  void processMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                 V0DatasLabeled const& V0Datas,
                 aod::McParticles const& mcparticles,
                 MyTracks const&)
  {

    // event selection
    if (collSelection && !collision.sel8()) {
      return;
    }
    hist.fill(HIST("hEventSelectionFlow"), 1.f);
    if (collSelection && (abs(collision.posZ()) >= 10.)) {
      return;
    }
    hist.fill(HIST("hEventSelectionFlow"), 2.f);

    for (auto& v0data : V0Datas) {

      if (v0data.has_mcParticle() && v0data.mcParticleId() > -1 && v0data.mcParticleId() <= mcparticles.size()) {
        auto MCv0 = v0data.mcParticle_as<aod::McParticles>();

        // Lambda
        if (MCv0.has_daughters() && MCv0.pdgCode() == LambdaPDG) {
          LOG(debug) << "V0 is a Lambda.";
          const auto& protonTrackIU = v0data.posTrack_as<MyTracks>();
          const auto& pionTrackIU = v0data.negTrack_as<MyTracks>();

          if (protonTrackIU.has_mcParticle() && pionTrackIU.has_mcParticle()) {

            const auto& MCproton = protonTrackIU.mcParticle_as<aod::McParticles>();
            const auto& MCpion = pionTrackIU.mcParticle_as<aod::McParticles>();

            if (MCproton.pdgCode() == ProtonPDG && MCpion.pdgCode() == NegPionPDG) {

              // lambda
              massLambda = v0data.mLambda();
              radiusLambda = v0data.v0radius();
              ptLambda = v0data.pt();
              /// daughter momenta at Lambda vertex
              // proton
              momProtonRec[0] = v0data.pxpos();
              momProtonRec[1] = v0data.pypos();
              momProtonRec[2] = v0data.pzpos();
              momProtonRecErr[0] = sqrt(v0data.covMatPosDau()[9]);
              momProtonRecErr[1] = sqrt(v0data.covMatPosDau()[14]);
              momProtonRecErr[2] = sqrt(v0data.covMatPosDau()[20]);
              momProtonGen[0] = MCproton.px();
              momProtonGen[1] = MCproton.py();
              momProtonGen[2] = MCproton.pz();
              // pion
              momPionRec[0] = v0data.pxneg();
              momPionRec[1] = v0data.pyneg();
              momPionRec[2] = v0data.pzneg();
              momPionRecErr[0] = sqrt(v0data.covMatNegDau()[9]);
              momPionRecErr[1] = sqrt(v0data.covMatNegDau()[14]);
              momPionRecErr[2] = sqrt(v0data.covMatNegDau()[20]);
              momPionGen[0] = MCpion.px();
              momPionGen[1] = MCpion.py();
              momPionGen[2] = MCpion.pz();

              // get daughter momenta at IU, charge, eta, nTPCclusters
              getTrackInfo(protonTrackIU, pionTrackIU);

              // fill table
              fillTable(collision);
            }
          }
        } // end Lambda

        // Anti-Lambda
        if (MCv0.pdgCode() == AntiLambdaPDG) {
          LOG(debug) << "V0 is an Anti-Lambda.";
          const auto& protonTrackIU = v0data.negTrack_as<MyTracks>();
          const auto& pionTrackIU = v0data.posTrack_as<MyTracks>();

          if (protonTrackIU.has_mcParticle() && pionTrackIU.has_mcParticle()) {

            const auto& MCproton = protonTrackIU.mcParticle_as<aod::McParticles>();
            const auto& MCpion = pionTrackIU.mcParticle_as<aod::McParticles>();

            if (MCproton.pdgCode() == AntiProtonPDG && MCpion.pdgCode() == PosPionPDG) {

              // lambda mass and radius
              massLambda = v0data.mAntiLambda();
              radiusLambda = v0data.v0radius();
              ptLambda = v0data.pt();
              /// daughter momenta at Lambda vertex
              // proton
              momProtonRec[0] = v0data.pxneg();
              momProtonRec[1] = v0data.pyneg();
              momProtonRec[2] = v0data.pzneg();
              momProtonRecErr[0] = sqrt(v0data.covMatNegDau()[9]);
              momProtonRecErr[1] = sqrt(v0data.covMatNegDau()[14]);
              momProtonRecErr[2] = sqrt(v0data.covMatNegDau()[20]);
              momProtonGen[0] = MCproton.px();
              momProtonGen[1] = MCproton.py();
              momProtonGen[2] = MCproton.pz();
              // pion
              momPionRec[0] = v0data.pxpos();
              momPionRec[1] = v0data.pypos();
              momPionRec[2] = v0data.pzpos();
              momPionRecErr[0] = sqrt(v0data.covMatPosDau()[9]);
              momPionRecErr[1] = sqrt(v0data.covMatPosDau()[14]);
              momPionRecErr[2] = sqrt(v0data.covMatPosDau()[20]);
              momPionGen[0] = MCpion.px();
              momPionGen[1] = MCpion.py();
              momPionGen[2] = MCpion.pz();

              // get daughter momenta at IU, charge, eta, nTPCclusters
              getTrackInfo(protonTrackIU, pionTrackIU);

              // fill table
              fillTable(collision);
            }
          }
        } // end Anti-Lambda
      }   // end MC
    }     // end V0 loop
  }
  PROCESS_SWITCH(qaLamMomResolution, processMC, "Process MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qaLamMomResolution>(cfgc)};
}
