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

/// \file HFXiOCCandidateSelector.cxx
/// \brief Omegac → Xi Pi selection task for Xi
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
using namespace o2::aod::cascdata;
using namespace o2::aod::cascdataext;
using namespace o2::aod::hf_cand_casc;
using namespace o2::aod::hf_cand_omegac;
using namespace o2::aod::hf_sel_omegac;
using namespace o2::aod::hf_lambdaoc;
using namespace o2::aod::hf_xioc;


/// Struct for applying Omegac selection cuts
struct HFXiOCCandidateSelector {
  Produces<aod::HFSelXiOC> hfSelXiOC;

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
  Configurable<double> cutinvmassxi{"cutinvmassxi", 0.05, "Invariant mass cut for xi (tot)"};

  //histograms
  OutputObj<TH1F> hSelections{TH1F("hSelections", "Sequential selections survivors;sequential selections;entries", 5,0.,5.)};
  OutputObj<TH1F> hInvMassXi{TH1F("hInvMassXi", "Xi invariant mass with invmass cut;inv mass;entries", 200, 1.27, 1.38)};
  OutputObj<TH1F> hInvMassXiCut{TH1F("hInvMassXiCut", "Xi invariant mass complete spectrum;inv mass;entries", 1000, 1.2, 1.4)};

  OutputObj<TH1F> hCheck1{TH1F("hCheck1", "Valies statuspidlambda;status pid lambda;entries", 5,-0.5,4.5)};
  OutputObj<TH1F> hCheck2{TH1F("hCheck2", "Valies statuspidantilambda;status pid antilambda;entries", 5,-0.5,4.5)};
  OutputObj<TH1F> hCheck3{TH1F("hCheck3", "Check inv mass casc initial;invmass;entries", 500,1.2,1.4)};
  OutputObj<TH1F> hCheck4{TH1F("hCheck4", "Values StatusInvMassXi;status;entries", 5,-0.5,4.5)};
  OutputObj<TH1F> hCheck5{TH1F("hCheck5", "Values StatusPidXi;status;entries", 5,-0.5,4.5)};

  //using FullV0Info = soa::Join<aod::V0Datas, aod::HFSelLambdaOC>;
  //aod::V0sLinked const&
  //aod::HFSelLambdaOC const&,

  void process(aod::CascDataExt const& candidates, aod::V0sLinked const&, aod::BigTracksPIDExtended const&)
    
  {
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(d_pidTPCMinpT, d_pidTPCMaxpT);
    selectorPion.setRangeNSigmaTPC(-d_nSigmaTPC, d_nSigmaTPC);
    selectorPion.setRangeNSigmaTPCCondTOF(-d_nSigmaTPCCombined, d_nSigmaTPCCombined);
    selectorPion.setRangePtTOF(d_pidTOFMinpT, d_pidTOFMaxpT);
    selectorPion.setRangeNSigmaTOF(-d_nSigmaTOF, d_nSigmaTOF);
    selectorPion.setRangeNSigmaTOFCondTPC(-d_nSigmaTOFCombined, d_nSigmaTOFCombined);

    // looping over lambda candidates
    for (auto& candidate : candidates) {

      auto trackPion = candidate.bachelor_as<aod::BigTracksPIDExtended>(); //pi <- Xi

      //-------------------- track-level PID selection --------------------
      int pidPion = selectorPion.getStatusTrackPIDTPC(trackPion);
      int StatusPidXi = 0;
      if(pidPion == TrackSelectorPID::Status::PIDAccepted){
        StatusPidXi = 1;
      }

      hCheck5->Fill(StatusPidXi);

      //-------------------- invariant mass cuts --------------------
      int StatusInvMassXi = 0;

      double invmassxi = -1; 
      if(StatusPidXi == 1){
        invmassxi = candidate.mXi();
      } 

      if (std::abs(invmassxi-1.32171) < (cutinvmassxi)){
        StatusInvMassXi = 1;
      }
      hCheck3->Fill(candidate.mXi());
      hCheck4->Fill(StatusInvMassXi);

     //filling
     hfSelXiOC(StatusPidXi, StatusInvMassXi);

     if(StatusPidXi == 1){
        hSelections->Fill(0.5);
        if(StatusInvMassXi == 1){
          hSelections->Fill(1.5);
          hInvMassXiCut->Fill(invmassxi);
        }
        hInvMassXi->Fill(invmassxi);
      }
      
    }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFXiOCCandidateSelector>(cfgc, TaskName{"hf-xioc-candidate-selector"})};
}
