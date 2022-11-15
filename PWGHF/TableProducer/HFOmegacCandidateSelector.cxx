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

/// \file HFOmegacCandidateSelector.cxx
/// \brief Omegac → Xi Pi selection task
///
/// \author Fedrica Zanone <federica.zanone@cern.ch>, Heidelberg University & GSI

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_casc;
using namespace o2::aod::hf_cand_omegac;
using namespace o2::aod::hf_sel_omegac;


/// Struct for applying Omegac selection cuts
struct HFOmegacCandidateSelector {
  Produces<aod::HFSelOmegacCandidate> hfSelOmegacCandidate;
  
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
  // topological cuts (see PWGHF/Core/HFSelectorCuts.h))
  //Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_d0_topik::pTBins_v}, "pT bin limits"};
  //Configurable<LabeledArray<double>> cuts{"D0_to_pi_K_cuts", {hf_cuts_d0_topik::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "D0 candidate selection per pT bin"};
  //invariant mass cuts
  Configurable<double> cutinvmasslambda{"cutinvmasslambda", 0.0025, "Invariant mass cut for lambda (sigma)"};
  Configurable<double> cutinvmasscascade{"cutinvmasscascade", 0.0025, "Invariant mass cut for cascade (sigma)"};
  Configurable<int> nsigmainvmasscut{"nsigmainvmasscut", 4, "Number of sigma for invariant mass cut"};
  Configurable<double> rangeinvmassomegac{"rangeinvmassomegac", 0.3, "Invariant mass range for omegac (sigma)"};

  OutputObj<TH1F> hPtPrimaryPi{TH1F("hPtPrimaryPi", "p_T primary #pi;p_T (GeV/#it{c});entries", 500, 0, 20)};
  OutputObj<TH1F> hxVertexOmegac{TH1F("hxVertexOmegac", "x Omegac vertex;xVtx;entries", 500, -10, 10)};
  OutputObj<TH1F> hInvMassLambda{TH1F("hInvMassLambda", "Lambda invariant mass;inv mass;entries", 1000, 1.10, 1.14)};
  OutputObj<TH1F> hInvMassCascade{TH1F("hInvMassCascade", "Cascade invariant mass;inv mass;entries", 1000, 1.295, 1.345)};
  OutputObj<TH1F> hInvMassOmegac{TH1F("hInvMassOmegac", "Omegac invariant mass;inv mass;entries", 500, 2.2, 3.1)};
  OutputObj<TH1F> hCTauOmegac{TH1F("hCTauOmegac", "Omegac ctau;ctau;entries", 500, 0., 10.)};

  OutputObj<TH1F> hTest1{TH1F("hTest1", "Test status steps;status;entries", 12, 0., 12.)};
  OutputObj<TH1F> hTest2{TH1F("hTest2", "Test status consecutive;status;entries", 12, 0., 12.)};



  /*
  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& candidate)
  {
    auto candpT = candidate.pt();
    auto pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < d_pTCandMin || candpT >= d_pTCandMax) {
      return false;
    }
    // product of daughter impact parameters
    if (candidate.impactParameterProduct() > cuts->get(pTBin, "d0d0")) {
      return false;
    }
    // cosine of pointing angle
    if (candidate.cpa() < cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }
    // cosine of pointing angle XY
    if (candidate.cpaXY() < cuts->get(pTBin, "cos pointing angle xy")) {
      return false;
    }
    // normalised decay length in XY plane
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normalized decay length XY")) {
      return false;
    }
    // candidate DCA
    //if (candidate.chi2PCA() > cuts[pTBin][1]) return false;

    // decay exponentail law, with tau = beta*gamma*ctau
    // decay length > ctau retains (1-1/e)
    if (std::abs(candidate.impactParameterNormalised0()) < 0.5 || std::abs(candidate.impactParameterNormalised1()) < 0.5) {
      return false;
    }
    double decayLengthCut = std::min((candidate.p() * 0.0066) + 0.01, cuts->get(pTBin, "minimum decay length"));
    if (candidate.decayLength() * candidate.decayLength() < decayLengthCut * decayLengthCut) {
      return false;
    }
    if (candidate.decayLength() > cuts->get(pTBin, "decay length")) {
      return false;
    }
    if (candidate.decayLengthXY() > cuts->get(pTBin, "decay length XY")) {
      return false;
    }
    if (candidate.decayLengthNormalised() * candidate.decayLengthNormalised() < 1.0) {
      //return false; // add back when getter fixed
    }
    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate is candidate
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \note trackPion = positive and trackKaon = negative for D0 selection and inverse for D0bar
  /// \return true if candidate passes all cuts for the given Conjugate
  template <typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackPion, const T2& trackKaon)
  {
    auto candpT = candidate.pt();
    auto pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // invariant-mass cut
    if (trackPion.sign() > 0) {
      if (std::abs(InvMassD0(candidate) - RecoDecay::getMassPDG(pdg::Code::kD0)) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(InvMassD0bar(candidate) - RecoDecay::getMassPDG(pdg::Code::kD0)) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    // cut on daughter pT
    if (trackPion.pt() < cuts->get(pTBin, "pT Pi") || trackKaon.pt() < cuts->get(pTBin, "pT K")) {
      return false;
    }

    // cut on daughter DCA - need to add secondary vertex constraint here
    if (std::abs(trackPion.dcaXY()) > cuts->get(pTBin, "d0pi") || std::abs(trackKaon.dcaXY()) > cuts->get(pTBin, "d0K")) {
      return false;
    }

    // cut on cos(theta*)
    if (trackPion.sign() > 0) {
      if (std::abs(CosThetaStarD0(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    } else {
      if (std::abs(CosThetaStarD0bar(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    }

    return true;
  }
  */
  
  void process(aod::HfCandOmegacBase const& candidates, aod::BigTracksPIDExtended const&)
  {
    TrackSelectorPID selectorPionFromOme(kPiPlus);
    selectorPionFromOme.setRangePtTPC(d_pidTPCMinpT, d_pidTPCMaxpT);
    selectorPionFromOme.setRangeNSigmaTPC(-d_nSigmaTPC, d_nSigmaTPC);
    selectorPionFromOme.setRangeNSigmaTPCCondTOF(-d_nSigmaTPCCombined, d_nSigmaTPCCombined);
    selectorPionFromOme.setRangePtTOF(d_pidTOFMinpT, d_pidTOFMaxpT);
    selectorPionFromOme.setRangeNSigmaTOF(-d_nSigmaTOF, d_nSigmaTOF);
    selectorPionFromOme.setRangeNSigmaTOFCondTPC(-d_nSigmaTOFCombined, d_nSigmaTOFCombined);

    TrackSelectorPID selectorPionFromCasc(kPiPlus);
    selectorPionFromCasc.setRangePtTPC(d_pidTPCMinpT, d_pidTPCMaxpT);
    selectorPionFromCasc.setRangeNSigmaTPC(-d_nSigmaTPC, d_nSigmaTPC);
    selectorPionFromCasc.setRangeNSigmaTPCCondTOF(-d_nSigmaTPCCombined, d_nSigmaTPCCombined);
    selectorPionFromCasc.setRangePtTOF(d_pidTOFMinpT, d_pidTOFMaxpT);
    selectorPionFromCasc.setRangeNSigmaTOF(-d_nSigmaTOF, d_nSigmaTOF);
    selectorPionFromCasc.setRangeNSigmaTOFCondTPC(-d_nSigmaTOFCombined, d_nSigmaTOFCombined);

    TrackSelectorPID selectorPionFromV0(kPiPlus);
    selectorPionFromV0.setRangePtTPC(d_pidTPCMinpT, d_pidTPCMaxpT);
    selectorPionFromV0.setRangeNSigmaTPC(-d_nSigmaTPC, d_nSigmaTPC);
    selectorPionFromV0.setRangeNSigmaTPCCondTOF(-d_nSigmaTPCCombined, d_nSigmaTPCCombined);
    selectorPionFromV0.setRangePtTOF(d_pidTOFMinpT, d_pidTOFMaxpT);
    selectorPionFromV0.setRangeNSigmaTOF(-d_nSigmaTOF, d_nSigmaTOF);
    selectorPionFromV0.setRangeNSigmaTOFCondTPC(-d_nSigmaTOFCombined, d_nSigmaTOFCombined);

    TrackSelectorPID selectorProton(kProton);
    selectorProton.setRangePtTPC(d_pidTPCMinpT, d_pidTPCMaxpT);
    selectorProton.setRangeNSigmaTPC(-d_nSigmaTPC, d_nSigmaTPC);
    selectorProton.setRangeNSigmaTPCCondTOF(-d_nSigmaTPCCombined, d_nSigmaTPCCombined);
    selectorProton.setRangePtTOF(d_pidTOFMinpT, d_pidTOFMaxpT);
    selectorProton.setRangeNSigmaTOF(-d_nSigmaTOF, d_nSigmaTOF);
    selectorProton.setRangeNSigmaTOFCondTPC(-d_nSigmaTOFCombined, d_nSigmaTOFCombined);

    // looping over omegac candidates
    for (auto& candidate : candidates) {

      // final selection flag: -1 - rejected, 1 - accepted

   

      /*if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        continue;
      }
      statusHFFlag = 1; */

      auto trackV0PosDau = candidate.posTrack_as<aod::BigTracksPIDExtended>(); // positive V0 daughter
      auto trackV0NegDau = candidate.negTrack_as<aod::BigTracksPIDExtended>(); // negative V0 daughter
      auto trackPiFromCasc = candidate.index2_as<aod::BigTracksPIDExtended>(); // pion <- cascade
      auto trackPiFromOmeg =  candidate.index1_as<aod::BigTracksPIDExtended>(); //pion <- omegac

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      int signdecay = candidate.signdecay();

      if(signdecay < 0){
        trackPiFromLam = trackV0PosDau;
        trackPrFromLam = trackV0NegDau;
      } else if (signdecay == 0){
        continue;
      }

      /*
      if (!daughterSelection(trackPos) || !daughterSelection(trackNeg)) {
        hfSelD0Candidate(statusD0, statusD0bar);
        continue;
      }
      */

      // conjugate-independent topological selection
      /*if (!selectionTopol(candidate)) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        continue;
      }
      statusTopol = 1; */

      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)

      // conjugate-dependent topological selection for D0
      //bool topolD0 = selectionTopolConjugate(candidate, trackPos, trackNeg);
      // conjugate-dependent topological selection for D0bar
      //bool topolD0bar = selectionTopolConjugate(candidate, trackNeg, trackPos);

      /*if (!topolD0 && !topolD0bar) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        continue;
      }
      statusCand = 1;*/

      // track-level PID selection
      int pidProton = selectorProton.getStatusTrackPIDTPC(trackPrFromLam);
      int pidPiFromLam = selectorPionFromV0.getStatusTrackPIDTPC(trackPiFromLam);
      int pidPiFromCasc = selectorPionFromCasc.getStatusTrackPIDTPC(trackPiFromCasc);
      int pidPiFromOme = selectorPionFromOme.getStatusTrackPIDTPC(trackPiFromOmeg);

      int statuspidLambda = -1;
      int statuspidCascade = -1;
      int statuspidOmegac = -1;

      if(pidProton == TrackSelectorPID::Status::PIDAccepted && pidPiFromLam == TrackSelectorPID::Status::PIDAccepted){
        statuspidLambda = 1;
        hTest2->Fill(0.5);
      }

      if(pidProton == TrackSelectorPID::Status::PIDAccepted && pidPiFromLam == TrackSelectorPID::Status::PIDAccepted && pidPiFromCasc == TrackSelectorPID::Status::PIDAccepted){
         statuspidCascade = 1;
         hTest2->Fill(1.5);
      }

      if(pidProton == TrackSelectorPID::Status::PIDAccepted && pidPiFromLam == TrackSelectorPID::Status::PIDAccepted && pidPiFromCasc == TrackSelectorPID::Status::PIDAccepted && pidPiFromOme == TrackSelectorPID::Status::PIDAccepted){
        statuspidOmegac = 1;
        hTest2->Fill(2.5);
      }

      // invariant mass cuts
      int statusinvmassLambda = -1;
      int statusinvmassCascade = -1;
      int statusinvmassOmegac = -1;

      double invmasslambda = 0; 
      if(signdecay > 0){
        invmasslambda = candidate.invmasslambda();
      } else if(signdecay < 0){
        invmasslambda = candidate.invmassantilambda();
      }
      double invmasscascade = candidate.invmasscascade();
      double invmassomegac = candidate.invmassomegac();

      if (abs(invmasslambda-1.11568) < (nsigmainvmasscut*cutinvmasslambda)){
        statusinvmassLambda = 1;
        if(statuspidLambda ==1 && statuspidCascade == 1 && statuspidOmegac == 1){
          hTest2->Fill(3.5);
        }
      }

      if (abs(invmasscascade-1.32171) < (nsigmainvmasscut*cutinvmasscascade)){
        statusinvmassCascade = 1;
        if(statuspidLambda ==1 && statuspidCascade == 1 && statuspidOmegac == 1 && statusinvmassLambda == 1){
          hTest2->Fill(4.5);
        }
      }

      if (abs(invmassomegac-2.6952) < rangeinvmassomegac){
        statusinvmassOmegac = 1;
        if(statuspidLambda ==1 && statuspidCascade == 1 && statuspidOmegac == 1 && statusinvmassLambda == 1 && statusinvmassCascade == 1){
          hTest2->Fill(5.5);
        };
      } 

      hfSelOmegacCandidate(statuspidLambda, statuspidCascade, statuspidOmegac, statusinvmassLambda, statusinvmassCascade, statusinvmassOmegac);
      
      if(statuspidLambda == -1){
        hTest1->Fill(0.5);
      }
      if(statuspidLambda == 1){
        hTest1->Fill(1.5);
      }
      if(statuspidCascade == -1){
        hTest1->Fill(2.5);
      }
      if(statuspidCascade == 1){
        hTest1->Fill(3.5);
      }
      if(statuspidOmegac == -1){
        hTest1->Fill(4.5);
      }
      if(statuspidOmegac == 1){
        hTest1->Fill(5.5);
      }
      if(statusinvmassLambda == -1){
        hTest1->Fill(6.5);
      }
      if(statusinvmassLambda == 1){
        hTest1->Fill(7.5);
      }
      if(statusinvmassCascade == -1){
        hTest1->Fill(8.5);
      }
      if(statusinvmassCascade == 1){
        hTest1->Fill(9.5);
      }
      if(statusinvmassOmegac == -1){
        hTest1->Fill(10.5);
      }
      if(statusinvmassOmegac == 1){
        hTest1->Fill(11.5);
      }



      double ptprimarypi = sqrt((candidate.pxprimarypiatprod()*candidate.pxprimarypiatprod())+(candidate.pyprimarypiatprod()*candidate.pyprimarypiatprod()));
      if(statuspidLambda == 1 && statuspidCascade== 1 && statuspidOmegac == 1 && statusinvmassLambda == 1 && statusinvmassCascade == 1 && statusinvmassOmegac == 1){
      //if(statuspidLambda == 1 && statuspidCascade== 1 && statuspidOmegac == 1  && statusinvmassOmegac == 1){
        hPtPrimaryPi->Fill(ptprimarypi);
        hxVertexOmegac->Fill(candidate.xdecayvtxomegac());
        hInvMassLambda->Fill(invmasslambda);
        hInvMassCascade->Fill(invmasscascade);
        hInvMassOmegac->Fill(invmassomegac);
        hCTauOmegac->Fill(candidate.ctauomegac());
      }
    }
  }
};






WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFOmegacCandidateSelector>(cfgc, TaskName{"hf-omegac-candidate-selector"})};
}
