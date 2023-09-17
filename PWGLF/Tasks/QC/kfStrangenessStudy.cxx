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
//
/// \brief this task allows for the direct one-to-one comparison of
//         cascades computed with standard DCAFitter methods and the KFparticle
//         package. It is meant for the purposes of larger-scale QA of KF reco.

/// \brief cascadebuilder.cxx task needs to be added to the workflow. Flag createCascCovMats needs to be enabled!

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/kfStrangenessStudy.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// allows for candidate-by-candidate comparison using Cascade to []CascData link table
using CascadesCrossLinked = soa::Join<aod::Cascades, aod::CascDataLink, aod::KFCascDataLink>;
using CascadesCrossLinkedLabeled = soa::Join<aod::Cascades, aod::CascDataLink, aod::KFCascDataLink, aod::McCascLabels>;

struct kfStrangenessStudy {

  Produces<aod::CascCand> rowCasc;
  Produces<aod::CascCandMC> rowCascMC;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// Filters
  Filter collisionFilter = (aod::evsel::sel8 == true);

  /// Cascade data
  int isDCAfitter = 0, isKF = 0;
  float pt = -1.0f, ptKF = -1.0f;
  float massXi = -1.0f, massXiKF = -1.0f;
  float massOmega = -1.0f, massOmegaKF = -1.0f;
  float massLambda = -1.0f, massLambdaKF = -1.0f;
  float cascRad = -1.0f, cascRadKF = -1.0f;
  float vtxX = -1.0f, vtxY = -1.0f, vtxZ = -1.0f, vtxXKF = -1.0f, vtxYKF = -1.0f, vtxZKF = -1.0f;
  float vtxXErr = -1.0f, vtxYErr = -1.0f, vtxZErr = -1.0f, vtxXErrKF = -1.0f, vtxYErrKF = -1.0f, vtxZErrKF = -1.0f;
  float dcaXYCascToPV = -1.0f, dcaXYCascToPVKF = -1.0f;
  float dcaZCascToPV = -1.0f, dcaZCascToPVKF = -1.0f;
  float dcaCascDaughters = -1.0f, dcaCascDaughtersKF = -1.0f;
  float dcaV0Daughters = -1.0f, dcaV0DaughtersKF = -1.0f;
  float dcaPosToPV = -1.0f, dcaPosToPVKF = -1.0f;
  float dcaNegToPV = -1.0f, dcaNegToPVKF = -1.0f;
  float dcaBachToPV = -1.0f, dcaBachToPVKF = -1.0f;
  float cascPointingAngle = -1.0f, cascPointingAngleKF = -1.0f;
  float V0Rad = -1.0f, V0RadKF = -1.0f;
  int charge = 0;
  float cascChi2geoKF = -1.f;

  /// Additional cascade MC data
  int isTrueCasc = 0;
  float ptRec = -1.0f, ptRecKF = -1.0f, ptGen = -1.0f;
  float vtxXrec = -1.0f, vtxYrec = -1.0f, vtxZrec = -1.0f, vtxXrecKF = -1.0f, vtxYrecKF = -1.0f, vtxZrecKF = -1.0f;
  float vtxXrecErr = -1.0f, vtxYrecErr = -1.0f, vtxZrecErr = -1.0f, vtxXrecErrKF = -1.0f, vtxYrecErrKF = -1.0f, vtxZrecErrKF = -1.0f;
  float vtxXgen = -1.0f, vtxYgen = -1.0f, vtxZgen = -1.0f;

  void init(InitContext const&)
  {
    /// QA histos
    histos.add("hChargeCounter", "hChargeCounter", kTH1F, {{3, -1.5f, 1.5f}});
    auto hEventSelectionFlow = histos.add<TH1>("hEventSelectionFlow", "Event selection flow", kTH1F, {{2, 0.5f, 2.5f}});
    hEventSelectionFlow->GetXaxis()->SetBinLabel(hEventSelectionFlow->FindBin(1), "Sel8");
    hEventSelectionFlow->GetXaxis()->SetBinLabel(hEventSelectionFlow->FindBin(2), "|Vtx_{z}|<10cm");
  }

  template <typename TCollision, typename TCascTable, typename TCascDatas, typename TKFCascDatas>
  void getCascDatas(TCollision const& collision, TCascTable const& cascade, TCascDatas const& cascdatas, TKFCascDatas const& kfcascdatas)
  {
    if (cascade.has_cascData()) {
      // check aod::Cascades -> aod::CascData link
      // if present: this candidate was accepted by default DCAfitter building
      isDCAfitter = 1;
      auto cascdata = cascade.template cascData_as<TCascDatas>();
      pt = cascdata.pt();
      massLambda = cascdata.mLambda();
      massXi = cascdata.mXi();
      massOmega = cascdata.mOmega();
      dcaXYCascToPV = cascdata.dcaXYCascToPV();
      dcaZCascToPV = cascdata.dcaZCascToPV();
      dcaCascDaughters = cascdata.dcacascdaughters();
      dcaV0Daughters = cascdata.dcaV0daughters();
      dcaPosToPV = cascdata.dcapostopv();
      dcaNegToPV = cascdata.dcanegtopv();
      dcaBachToPV = cascdata.dcabachtopv();
      cascPointingAngle = TMath::ACos(cascdata.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      cascRad = cascdata.cascradius();
      vtxX = cascdata.x();
      vtxY = cascdata.y();
      vtxZ = cascdata.z();
      vtxXErr = sqrt(cascdata.positionCovMat()[0]);
      vtxYErr = sqrt(cascdata.positionCovMat()[2]);
      vtxZErr = sqrt(cascdata.positionCovMat()[5]);
      V0Rad = cascdata.v0radius();
    }
    if (cascade.has_kfCascData()) {
      // check aod::Cascades -> aod::KFCascData link
      // if present: this candidate was accepted by KF building
      isKF = 1;
      auto cascdata = cascade.template kfCascData_as<TKFCascDatas>();
      ptKF = cascdata.pt();
      massLambdaKF = cascdata.mLambda();
      massXiKF = cascdata.mXi();
      massOmegaKF = cascdata.mOmega();
      dcaXYCascToPVKF = cascdata.dcaXYCascToPV();
      dcaZCascToPVKF = cascdata.dcaZCascToPV();
      dcaCascDaughtersKF = cascdata.dcacascdaughters();
      dcaV0DaughtersKF = cascdata.dcaV0daughters();
      dcaPosToPVKF = cascdata.dcapostopv();
      dcaNegToPVKF = cascdata.dcanegtopv();
      dcaBachToPVKF = cascdata.dcabachtopv();
      cascPointingAngleKF = TMath::ACos(cascdata.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      cascRadKF = cascdata.cascradius();
      vtxXKF = cascdata.x();
      vtxYKF = cascdata.y();
      vtxZKF = cascdata.z();
      vtxXErrKF = sqrt(cascdata.kfTrackCovMat()[0]);
      vtxYErrKF = sqrt(cascdata.kfTrackCovMat()[2]);
      vtxZErrKF = sqrt(cascdata.kfTrackCovMat()[5]);
      V0RadKF = cascdata.v0radius();
      cascChi2geoKF = cascdata.kfCascadeChi2();
    }
  }

  template <typename TCollision, typename TCascTable>
  void fillCascDataTable(TCollision const& collision, TCascTable const& cascade)
  {
    rowCasc(collision.globalIndex(),
          pt, ptKF,
          massXi, massXiKF,
          massOmega, massOmegaKF,
          cascRad, cascRadKF,
          vtxX, vtxY, vtxZ, vtxXErr, vtxYErr, vtxZErr,
          vtxXKF, vtxXKF, vtxXKF, vtxXErrKF, vtxYErrKF, vtxZErrKF,
          dcaXYCascToPV, dcaXYCascToPVKF,
          dcaZCascToPV, dcaZCascToPVKF,
          dcaCascDaughters, dcaCascDaughtersKF,
          cascPointingAngle, cascPointingAngleKF,
          charge, 
          massLambda, massLambdaKF,
          V0Rad, V0RadKF,
          dcaV0Daughters, dcaV0DaughtersKF,
          dcaPosToPV, dcaPosToPVKF,
          dcaNegToPV, dcaNegToPVKF,
          dcaBachToPV, dcaBachToPVKF,
          isDCAfitter, isKF);
  }

  template <typename TCollision, typename TCascTable>
  void fillCascMCTable(TCollision const& collision, TCascTable const& cascade)
  {
    rowCascMC(collision.globalIndex(),
            ptRec, ptRecKF, ptGen,
            massXi, massXiKF,
            massOmega, massOmegaKF,
            cascRad, cascRadKF,
            vtxXrec, vtxYrec, vtxZrec, vtxXrecErr, vtxYrecErr, vtxZrecErr,
            vtxXrecKF, vtxXrecKF, vtxXrecKF, vtxXrecErrKF, vtxYrecErrKF, vtxZrecErrKF,
            vtxXgen, vtxYgen, vtxZgen,
            dcaXYCascToPV, dcaXYCascToPVKF,
            dcaZCascToPV, dcaZCascToPVKF,
            dcaCascDaughters, dcaCascDaughtersKF,
            cascPointingAngle, cascPointingAngleKF,
            charge, 
            massLambda, massLambdaKF,
            V0Rad, V0RadKF,
            dcaV0Daughters, dcaV0DaughtersKF,
            dcaPosToPV, dcaPosToPVKF,
            dcaNegToPV, dcaNegToPVKF,
            dcaBachToPV, dcaBachToPVKF,
            isDCAfitter, isKF,
            isTrueCasc);
  }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, CascadesCrossLinked const& Cascades, soa::Join<aod::CascDatas, aod::CascCovs> const& CascDatas, soa::Join<aod::KFCascDatas, aod::KFCascCovs> const& KFCascDatas, aod::TracksIU const&)
  {
    /// Event selection
    histos.fill(HIST("hEventSelectionFlow"), 1.f);
    // select collisions in acceptance
    if (!(abs(collision.posZ()) < 10.)) return;
    histos.fill(HIST("hEventSelectionFlow"), 2.f);
   
   for (auto& cascade : Cascades) { // allows for cross-referencing everything

      // get charge from bachelor (unambiguous wrt to building)
      auto bachTrack = cascade.bachelor_as<aod::TracksIU>();
      if (bachTrack.sign() < 0) {
        charge = -1;
      } else {
        charge = +1;
      }
      histos.fill(HIST("hChargeCounter"), charge);

      // get cascade data and fill table
      getCascDatas(collision, cascade, CascDatas, KFCascDatas);
      if (cascade.has_cascData() || cascade.has_kfCascData()) {
        fillCascDataTable(collision, cascade);
      }
    } // end cascade loop
  } // end process
  PROCESS_SWITCH(kfStrangenessStudy, processData, "process data", true);

  void processMC(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, CascadesCrossLinkedLabeled const& Cascades, soa::Join<aod::CascDatas, aod::CascCovs> const& CascDatas, soa::Join<aod::KFCascDatas, aod::KFCascCovs> const& KFCascDatas,  aod::TracksIU const&, aod::McParticles const& particlesMC) 
  {
    /// Event selection
    histos.fill(HIST("hEventSelectionFlow"), 1.f);
    // select collisions in acceptance
    if (!(abs(collision.posZ()) < 10.)) return;
    histos.fill(HIST("hEventSelectionFlow"), 2.f);

    for (auto& cascade : Cascades) { // allows for cross-referencing everything

      // get charge from bachelor (unambiguous wrt to building)
      auto bachTrack = cascade.bachelor_as<aod::TracksIU>();
      if (bachTrack.sign() < 0) {
        charge = -1;
      } else {
        charge = +1;
      }
      histos.fill(HIST("hChargeCounter"), charge);

      // get cascade data
      getCascDatas(collision, cascade, CascDatas, KFCascDatas);

      // get cascade MC information
      if (cascade.has_cascData() || cascade.has_kfCascData()) {
        if (cascade.mcParticleId() >= -1 && cascade.mcParticleId() <= particlesMC.size()) {
          auto MCcascade = cascade.mcParticle_as<aod::McParticles>();
          ptGen = MCcascade.pt();
          vtxXgen = MCcascade.vx();
          vtxYgen = MCcascade.vy();
          vtxZgen = MCcascade.vz();
          if (abs(MCcascade.pdgCode()) == 3312) { // Xi
            isTrueCasc = 1;
          } else {
            isTrueCasc = 0;
          }
        }

        // fill cascade table
        fillCascMCTable(collision, cascade);

      }
    } // end cascade loop
  } // end process
  PROCESS_SWITCH(kfStrangenessStudy, processMC, "process MC", false);

}; 

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kfStrangenessStudy>(cfgc)};
}
