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

struct kfStrangenessStudy {

  Produces<aod::CascCand> rowCasc;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    histos.add("hEventCounter", "hEventCounter", kTH1F, {{1, 0.0f, 1.0f}});
    histos.add("hChargeCounter", "hChargeCounter", kTH1F, {{3, -1.5f, 1.5f}});
  }

  void process(aod::Collision const& collision, CascadesCrossLinked const& Cascades, aod::CascDatas const&, aod::KFCascDatas const&, aod::TracksIU const&)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& cascade : Cascades) { // allows for cross-referencing everything
      int isDCAfitter = 0, isKF = 0;
      float pt = 0.0f, ptKF = 0.0f;
      float massXi = 0.0f, massXiKF = 0.0f;
      float massOmega = 0.0f, massOmegaKF = 0.0f;
      float massLambda = 0.0f, massLambdaKF = 0.0f;
      float dcaXYCascToPV = 0.0f, dcaXYCascToPVKF = 0.0f;
      float dcaZCascToPV = 0.0f, dcaZCascToPVKF = 0.0f;
      float dcaCascDaughters = 0.0f, dcaCascDaughtersKF = 0.0f;
      // float dcaV0ToPV = 0.0f, dcaV0ToPVKF = 0.0f;
      float dcaV0Daughters = 0.0f, dcaV0DaughtersKF = 0.0f;
      float dcaPosToPV = 0.0f, dcaPosToPVKF = 0.0f;
      float dcaNegToPV = 0.0f, dcaNegToPVKF = 0.0f;
      float dcaBachToPV = 0.0f, dcaBachToPVKF = 0.0f;
      float cascPointingAngle = -1.0f, cascPointingAngleKF = -1.0f;
      float cascRad = 0.0f, cascRadKF = 0.0f;
      float V0Rad = 0.0f, V0RadKF = 0.0f;
      int charge = 0;

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
        // dcaV0ToPV = cascdata.dcav0topv();
        dcaV0Daughters = cascdata.dcaV0daughters();
        dcaPosToPV = cascdata.dcapostopv();
        dcaNegToPV = cascdata.dcanegtopv();
        dcaBachToPV = cascdata.dcabachtopv();
        cascPointingAngle = TMath::ACos(cascdata.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        cascRad = cascdata.cascradius();
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
        // dcaV0ToPVKF = cascdata.dcav0topv();
        dcaV0DaughtersKF = cascdata.dcaV0daughters();
        dcaPosToPVKF = cascdata.dcapostopv();
        dcaNegToPVKF = cascdata.dcanegtopv();
        dcaBachToPVKF = cascdata.dcabachtopv();
        cascPointingAngleKF = TMath::ACos(cascdata.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        cascRadKF = cascdata.cascradius();
        V0RadKF = cascdata.v0radius();
      }

      /// write cascade tree
      rowCasc(collision.globalIndex(),
              pt, ptKF,
              massXi, massXiKF,
              massOmega, massOmegaKF,
              massLambda, massLambdaKF,
              dcaXYCascToPV, dcaXYCascToPVKF,
              dcaZCascToPV, dcaZCascToPVKF,
              dcaCascDaughters, dcaCascDaughtersKF,
              //dcaV0ToPV, dcaV0ToPVKF,
              dcaV0Daughters, dcaV0DaughtersKF,
              dcaPosToPV, dcaPosToPVKF,
              dcaNegToPV, dcaNegToPVKF,
              dcaBachToPV, dcaBachToPVKF,
              cascPointingAngle, cascPointingAngleKF,
              cascRad, cascRadKF,
              V0Rad, V0RadKF,
              charge, 
              isDCAfitter, isKF);
    } // end cascade loop
  } // end process
}; 

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kfStrangenessStudy>(cfgc)};
}
