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

/// \file HFLambdaOCCandidateSelector.cxx
/// \brief Omegac → Xi Pi selection task for lambda
///
/// \author Fedrica Zanone <federica.zanone@cern.ch>, Heidelberg University & GSI

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_casc;
using namespace o2::aod::hf_cand_omegac;
using namespace o2::aod::hf_sel_omegac;
using namespace o2::aod::hf_lambdaoc;


/// Struct for applying Omegac selection cuts
struct HFLambdaOCCandidateSelector {
  Produces<aod::HFSelLambdaOC> hfSelLambdaOC;

  Configurable<double> d_pTCandMin{"d_pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> d_pTCandMax{"d_pTCandMax", 50., "Upper bound of candidate pT"};
  // TPC
  Configurable<double> d_pidTPCMinpT{"d_pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> d_pidTPCMaxpT{"d_pidTPCMaxpT", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> d_nSigmaTPC{"d_nSigmaTPC", 3., "Nsigma cut on TPC only"};
  Configurable<double> d_nSigmaTPCCombined{"d_nSigmaTPCCombined", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF
  Configurable<double> d_pidTOFMinpT{"d_pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> d_pidTOFMaxpT{"d_pidTOFMaxpT", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> d_nSigmaTOF{"d_nSigmaTOF", 3., "Nsigma cut on TOF only"};
  Configurable<double> d_nSigmaTOFCombined{"d_nSigmaTOFCombined", 5., "Nsigma cut on TOF combined with TPC"};
  //invariant mass cuts
  Configurable<double> cutinvmasslambda{"cutinvmasslambda", 0.05, "Invariant mass cut for lambda (tot)"};

  //histograms
  OutputObj<TH1F> hSelections{TH1F("hSelections", "Sequential selections survivors;sequential selections;entries", 5,0.,5.)};
  OutputObj<TH1F> hInvMassLambda{TH1F("hInvMassLambda", "Lambda invariant mass with invmass cut;inv mass;entries", 500, 1.06, 1.17)};
  OutputObj<TH1F> hInvMassLambdaCut{TH1F("hInvMassLambdaCut", "Lambda invariant mass complete spectrum;inv mass;entries", 1000, 1.0, 1.2)};

  OutputObj<TH1F> hCheck1{TH1F("hCheck1", "Valies statuspidlambda;status pid lambda;entries", 5,-0.5,4.5)};
  OutputObj<TH1F> hCheck2{TH1F("hCheck2", "Valies statuspidantilambda;status pid antilambda;entries", 5,-0.5,4.5)};


  
  void process(aod::V0Datas const& candidates, aod::V0sLinked const&, aod::BigTracksPIDExtended const&)
  {
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(d_pidTPCMinpT, d_pidTPCMaxpT);
    selectorPion.setRangeNSigmaTPC(-d_nSigmaTPC, d_nSigmaTPC);
    selectorPion.setRangeNSigmaTPCCondTOF(-d_nSigmaTPCCombined, d_nSigmaTPCCombined);
    selectorPion.setRangePtTOF(d_pidTOFMinpT, d_pidTOFMaxpT);
    selectorPion.setRangeNSigmaTOF(-d_nSigmaTOF, d_nSigmaTOF);
    selectorPion.setRangeNSigmaTOFCondTPC(-d_nSigmaTOFCombined, d_nSigmaTOFCombined);

    TrackSelectorPID selectorProton(kProton);
    selectorProton.setRangePtTPC(d_pidTPCMinpT, d_pidTPCMaxpT);
    selectorProton.setRangeNSigmaTPC(-d_nSigmaTPC, d_nSigmaTPC);
    selectorProton.setRangeNSigmaTPCCondTOF(-d_nSigmaTPCCombined, d_nSigmaTPCCombined);
    selectorProton.setRangePtTOF(d_pidTOFMinpT, d_pidTOFMaxpT);
    selectorProton.setRangeNSigmaTOF(-d_nSigmaTOF, d_nSigmaTOF);
    selectorProton.setRangeNSigmaTOFCondTPC(-d_nSigmaTOFCombined, d_nSigmaTOFCombined);

    // looping over lambda candidates
    for (auto& candidate : candidates) {

      auto trackV0PosDau = candidate.posTrack_as<aod::BigTracksPIDExtended>(); // positive V0 daughter
      auto trackV0NegDau = candidate.negTrack_as<aod::BigTracksPIDExtended>(); // negative V0 daughter

      //-------------------- track-level PID selection --------------------
      //decay of lambda
      int pidProton = selectorProton.getStatusTrackPIDTPC(trackV0PosDau);
      int pidPionMinus = selectorPion.getStatusTrackPIDTPC(trackV0NegDau);
      //decay of anti-lambda
      int pidAntiProton = selectorProton.getStatusTrackPIDTPC(trackV0NegDau);
      int pidPionPlus = selectorPion.getStatusTrackPIDTPC(trackV0PosDau);

      int StatusPidLambda = 0;
      //0 -> PID failed
      //1 -> PID succeeded for lambda decay
      int StatusPidAntiLambda = 0;
      //0 -> PID failed
      //1 -> PID succeeded for antilambda decay

      if(pidProton == TrackSelectorPID::Status::PIDAccepted && pidPionMinus == TrackSelectorPID::Status::PIDAccepted){
        StatusPidLambda = 1;
      }

      if(pidAntiProton == TrackSelectorPID::Status::PIDAccepted && pidPionPlus == TrackSelectorPID::Status::PIDAccepted){
         StatusPidAntiLambda = 1;
      }

      //-------------------- invariant mass cuts --------------------
      int StatusInvMassLambda = 0;

      double invmasslambda = -1; 
      if(StatusPidLambda == 1){
        invmasslambda = candidate.mLambda();
      } else if(StatusPidAntiLambda == 1){
        invmasslambda = candidate.mAntiLambda();
      } else if(StatusPidLambda == 1 && StatusPidAntiLambda){
        invmasslambda = -1;
      }

      if (std::abs(invmasslambda-1.11568) < (cutinvmasslambda)){
        StatusInvMassLambda = 1;
      }

      hfSelLambdaOC(candidate.globalIndex(), StatusPidLambda, StatusPidAntiLambda, StatusInvMassLambda);

      hCheck1->Fill(StatusPidLambda);
      hCheck2->Fill(StatusPidAntiLambda);

     if(!(StatusPidLambda == 1 && StatusPidAntiLambda == 1)){
     if(StatusPidLambda == 1 || StatusPidAntiLambda == 1){
        hSelections->Fill(0.5);
        if(StatusInvMassLambda == 1){
          hSelections->Fill(1.5);
          hInvMassLambdaCut->Fill(invmasslambda);
        }
        hInvMassLambda->Fill(invmasslambda);
      }
     }
      
    }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFLambdaOCCandidateSelector>(cfgc, TaskName{"hf-lambdaoc-candidate-selector"})};
}
