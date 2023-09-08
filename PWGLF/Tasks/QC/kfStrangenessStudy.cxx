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
using CascadesCrossLinked = soa::Join<aod::Cascades, aod::CascDataLink, aod::KFCascDataLink, o2::aod::CascCovs, aod::KFCascCovs>;
using CascadesCrossLinkedLabeled = soa::Join<aod::Cascades, aod::CascDataLink, aod::KFCascDataLink, o2::aod::CascCovs, aod::KFCascCovs, aod::McCascLabels>;

struct kfStrangenessStudy {

  Produces<aod::CascCand> rowCasc;
  Produces<aod::CascCandMC> rowCascMC;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// Filters
  Filter collisionFilter = (aod::evsel::sel8 == true);

  /// Cascade data
  int isDCAfitter = 0, isKF = 0;
  float pt = 0.0f, ptKF = 0.0f;
  float massXi = 0.0f, massXiKF = 0.0f;
  float massOmega = 0.0f, massOmegaKF = 0.0f;
  float massLambda = 0.0f, massLambdaKF = 0.0f;
  float cascRad = 0.0f, cascRadKF = 0.0f;
  float vtxX = 0.0f, vtxY = 0.0f, vtxZ = 0.0f, vtxXKF = 0.0f, vtxYKF = 0.0f, vtxZKF = 0.0f;
  float vtxXErr = 0.0f, vtxYErr = 0.0f, vtxZErr = 0.0f, vtxXErrKF = 0.0f, vtxYErrKF = 0.0f, vtxZErrKF = 0.0f;
  float dcaXYCascToPV = 0.0f, dcaXYCascToPVKF = 0.0f;
  float dcaZCascToPV = 0.0f, dcaZCascToPVKF = 0.0f;
  float dcaCascDaughters = 0.0f, dcaCascDaughtersKF = 0.0f;
  float dcaV0Daughters = 0.0f, dcaV0DaughtersKF = 0.0f;
  float dcaPosToPV = 0.0f, dcaPosToPVKF = 0.0f;
  float dcaNegToPV = 0.0f, dcaNegToPVKF = 0.0f;
  float dcaBachToPV = 0.0f, dcaBachToPVKF = 0.0f;
  float cascPointingAngle = -1.0f, cascPointingAngleKF = -1.0f;
  float V0Rad = 0.0f, V0RadKF = 0.0f;
  int charge = 0;

  /// Additional cascade MC data
  int isTrueCasc = 0;
  float ptRec = 0.0f, ptRecKF = 0.0f, ptGen = 0.0f;
  float vtxXrec = 0.0f, vtxYrec = 0.0f, vtxZrec = 0.0f, vtxXrecKF = 0.0f, vtxYrecKF = 0.0f, vtxZrecKF = 0.0f;
  float vtxXrecErr = 0.0f, vtxYrecErr = 0.0f, vtxZrecErr = 0.0f, vtxXrecErrKF = 0.0f, vtxYrecErrKF = 0.0f, vtxZrecErrKF = 0.0f;
  float vtxXgen = 0.0f, vtxYgen = 0.0f, vtxZgen = 0.0f;

  void init(InitContext const&)
  {
    /// QA histos
    histos.add("hChargeCounter", "hChargeCounter", kTH1F, {{3, -1.5f, 1.5f}});
    auto hEventSelectionFlow = histos.add<TH1>("hEventSelectionFlow", "Event selection flow", kTH1F, {{2, 0.5f, 2.5f}});
    hEventSelectionFlow->GetXaxis()->SetBinLabel(hEventSelectionFlow->FindBin(1), "Sel8");
    hEventSelectionFlow->GetXaxis()->SetBinLabel(hEventSelectionFlow->FindBin(2), "|Vtx_{z}|<10cm");
  }

  void processData(soa::Filtered<aod::Collision>::iterator const& collision, CascadesCrossLinked const& Cascades, aod::CascDatas const&, aod::KFCascDatas const&, aod::TracksIU const&)
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

      if (cascade.has_cascData()) {
        // check aod::Cascades -> aod::CascData link
        // if present: this candidate was accepted by default DCAfitter building
        isDCAfitter = 1;
        auto cascdata = cascade.cascData();
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
        vtxXErr = sqrt(cascade.positionCovMat()[0]);
        vtxYErr = sqrt(cascade.positionCovMat()[2]);
        vtxZErr = sqrt(cascade.positionCovMat()[5]);
        V0Rad = cascdata.v0radius();
      }
      if (cascade.has_kfCascData()) {
        // check aod::Cascades -> aod::KFCascData link
        // if present: this candidate was accepted by KF building
        isKF = 1;
        auto cascdata = cascade.kfCascData();
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
        vtxXErrKF = sqrt(cascade.kfTrackCovMat()[0]);
        vtxYErrKF = sqrt(cascade.kfTrackCovMat()[2]);
        vtxZErrKF = sqrt(cascade.kfTrackCovMat()[5]);
        V0RadKF = cascdata.v0radius();
      }

      /// write cascade tree
      if (cascade.has_kfCascData() || cascade.has_kfCascData()) {
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
    } // end cascade loop
  } // end process
  PROCESS_SWITCH(kfStrangenessStudy, processData, "process data", true);

  void processMC(soa::Filtered<aod::Collision>::iterator const& collision, CascadesCrossLinkedLabeled const& Cascades, aod::CascDatas const&, aod::KFCascDatas const&, aod::TracksIU const&, aod::McParticles const& particlesMC) 
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

      if (cascade.has_cascData()) {
        // check aod::Cascades -> aod::CascData link
        // if present: this candidate was accepted by default DCAfitter building
        isDCAfitter = 1;
        auto cascdata = cascade.cascData();
        ptRec = cascdata.pt();
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
        vtxXrec = cascdata.x();
        vtxYrec = cascdata.y();
        vtxZrec = cascdata.z();
        vtxXrecErr = sqrt(cascade.positionCovMat()[0]);
        vtxYrecErr = sqrt(cascade.positionCovMat()[2]);
        vtxZrecErr = sqrt(cascade.positionCovMat()[5]);
        V0Rad = cascdata.v0radius();

      }
      if (cascade.has_kfCascData()) {
        // check aod::Cascades -> aod::KFCascData link
        // if present: this candidate was accepted by KF building
        isKF = 1;
        auto cascdata = cascade.kfCascData();
        ptRecKF = cascdata.pt();
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
        vtxXrecKF = cascdata.x();
        vtxYrecKF = cascdata.y();
        vtxZrecKF = cascdata.z();
        vtxXrecErrKF = sqrt(cascade.kfTrackCovMat()[0]);
        vtxYrecErrKF = sqrt(cascade.kfTrackCovMat()[2]);
        vtxZrecErrKF = sqrt(cascade.kfTrackCovMat()[5]);
        V0RadKF = cascdata.v0radius();
      }

      if (cascade.has_kfCascData() || cascade.has_kfCascData()) {

        /// Get MC information
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

        /// write cascade tree
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
    } // end cascade loop
  } // end process
  PROCESS_SWITCH(kfStrangenessStudy, processMC, "process MC", false);

}; 

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kfStrangenessStudy>(cfgc)};
}
