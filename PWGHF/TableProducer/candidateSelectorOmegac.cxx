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
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Common/Core/TrackSelection.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_omegac;
using namespace o2::aod::hf_sel_omegac;

/// Struct for applying Omegac selection cuts
struct HFOmegacCandidateSelector {
  Produces<aod::HFSelOmegacCandidate> hfSelOmegacCandidate;

  // LF analysis selections
  // zPV -> can be already set in HFeventselection -> 10 cm
  // sel8 -> can be already set in HFeventselection -> true
  Configurable<double> CascRadius{"CascRadius", 0.5, "Min cascade radius"};
  Configurable<double> V0Radius{"V0Radius", 2, "Min V0 radius"};
  Configurable<double> V0CosPA{"V0CosPA", 0.95, "Cos PA V0"};
  Configurable<double> CascCosPA{"CascCosPA", 0.95, "Cos PA cascade"};
  Configurable<double> dcaCascDau{"dcaCascDau", 5.0, "Max DCA cascade daughters"};
  Configurable<double> dcaV0Dau{"dcaV0Dau", 5.0, "Max DCA V0 daughters"};
  Configurable<double> dcaOmegacDau{"dcaOmegacDau", 5.0, "Max DCA omegac daughters"};

  // limit charm baryon invariant mass spectrum
  Configurable<double> LowerLimitSpectrum{"LowerLimitSpectrum", 2.4, "Lower limit invariant mass spectrum charm baryon"};
  Configurable<double> UpperLimitSpectrum{"UpperLimitSpectrum", 3.0, "Upper limit invariant mass spectrum charm baryon"};

  // kinematic selections
  Configurable<double> EtaMax{"EtaMax", 0.8, "Max absolute value of eta"};
  Configurable<double> pTMinPiFromCasc{"pTMinPiFromCasc", 0.15, "Min pT pi <- casc"};
  Configurable<double> pTMinPiFromOme{"pTMinPiFromOme", 0.2, "Min pT pi <- omegac"};
  Configurable<double> PrimPiMindcaxyMin{"PrimPiMindcaxyMin", 0., "Min dcaxy primary pi track to PV"};
  Configurable<double> PrimPiMindcaxyMax{"PrimPiMindcaxyMax", 10., "Max dcaxy primary pi track to PV"};

  Configurable<double> d_pTCandMin{"d_pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> d_pTCandMax{"d_pTCandMax", 50., "Upper bound of candidate pT"};

  // PID options
  Configurable<bool> PIDTPCOnly{"PIDTPCOnly", true, "PID with TPC only"};
  Configurable<bool> PIDTPCTOFCombined{"PIDTPCTOFCombined", false, "PID with TPC & TOF"};

  // PID - TPC selections
  Configurable<double> d_pidTPCMinpT{"d_pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> d_pidTPCMaxpT{"d_pidTPCMaxpT", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> d_nSigmaTPC{"d_nSigmaTPC", 3., "Nsigma cut on TPC only"};
  Configurable<double> d_nSigmaTPCCombined{"d_nSigmaTPCCombined", 5., "Nsigma cut on TPC combined with TOF"};
  // PID - TOF selections
  Configurable<double> d_pidTOFMinpT{"d_pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> d_pidTOFMaxpT{"d_pidTOFMaxpT", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> d_nSigmaTOF{"d_nSigmaTOF", 3., "Nsigma cut on TOF only"};
  Configurable<double> d_nSigmaTOFCombined{"d_nSigmaTOFCombined", 5., "Nsigma cut on TOF combined with TPC"};

  // invariant mass cuts
  Configurable<double> cutinvmasslambda{"cutinvmasslambda", 0.0025, "Invariant mass cut for lambda (sigma)"};
  Configurable<double> cutinvmasscascade{"cutinvmasscascade", 0.0025, "Invariant mass cut for cascade (sigma)"};
  Configurable<int> nsigmainvmasscut{"nsigmainvmasscut", 4, "Number of sigma for invariant mass cut"};
  Configurable<double> rangeinvmassomegac{"rangeinvmassomegac", 0.3, "Invariant mass range for omegac (sigma)"};

  // detector clusters selections
  Configurable<int> tpcClusters{"tpcClusters", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> tpcCrossedRows{"tpcCrossedRows", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> tpcCrossedRowsOverFindableClustersRatio{"tpcCrossedRowsOverFindableClustersRatio", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<int> itsClusters{"itsClusters", 3, "Minimum number of ITS clusters requirement for pi <- Omegac"};
  Configurable<int> itsClustersInnBarr{"itsClustersInnBarr", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- Omegac"};
  // NOTA: MinNumberOfTPCCrossedRows=70 and MinTPCCrossedRows/FindableClustersRatio=0.8 already implemented in GlobalTrack selections (it's also required to have at least 1 hit in any SPD layer)

  OutputObj<TH1F> hPtPrimaryPi{TH1F("hPtPrimaryPi", "p_T primary #pi;p_T (GeV/#it{c});entries", 500, 0, 20)};
  OutputObj<TH1F> hxVertexOmegac{TH1F("hxVertexOmegac", "x Omegac vertex;xVtx;entries", 500, -10, 10)};
  OutputObj<TH1F> hInvMassOmegac{TH1F("hInvMassOmegac", "Omegac invariant mass;inv mass;entries", 500, 2.2, 3.1)};
  OutputObj<TH1F> hCTauOmegac{TH1F("hCTauOmegac", "Omegac ctau;ctau;entries", 500, 0., 10.)};
  OutputObj<TH1F> hInvMassOmegacNotFixed{TH1F("hInvMassOmegacNotFixed", "Omegac invariant mass (not fixed);inv mass;entries", 500, 2.2, 3.1)};

  OutputObj<TH1F> hTest1{TH1F("hTest1", "Test status steps;status;entries", 12, 0., 12.)};
  OutputObj<TH1F> hTest2{TH1F("hTest2", "Test status consecutive;status;entries", 12, 0., 12.)};

  using MyTrackInfo = soa::Join<aod::BigTracksPIDExtended, aod::TrackSelection>;

  void process(aod::HfCandOmegacBase const& candidates, MyTrackInfo const&)
  // aod::BigTracksPIDExtended const&)
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

      auto trackV0PosDau = candidate.posTrack_as<MyTrackInfo>();    // positive V0 daughter
      auto trackV0NegDau = candidate.negTrack_as<MyTrackInfo>();    // negative V0 daughter
      auto trackPiFromCasc = candidate.bachelor_as<MyTrackInfo>();  // pion <- cascade
      auto trackPiFromOmeg = candidate.primarypi_as<MyTrackInfo>(); // pion <- omegac

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      int signdecay = candidate.signdecay(); // sign of pi <- cascade

      if (signdecay > 0) {
        trackPiFromLam = trackV0PosDau;
        trackPrFromLam = trackV0NegDau;
      } else if (signdecay == 0) {
        continue;
      }

      // eta selection
      double etav0posdau = candidate.etav0posdau();
      double etav0negdau = candidate.etav0negdau();
      double etapifromcasc = candidate.etapifromcasc();
      double etapifromome = candidate.etapifromome();
      if (abs(etav0posdau) > EtaMax) {
        continue;
      }
      if (abs(etav0negdau) > EtaMax) {
        continue;
      }
      if (abs(etapifromcasc) > EtaMax) {
        continue;
      }
      if (abs(etapifromome) > EtaMax) {
        continue;
      }

      // minimum radius cut (LFcut)
      if (RecoDecay::sqrtSumOfSquares(candidate.xdecayvtxcascade(), candidate.ydecayvtxcascade()) < CascRadius) {
        continue;
      }
      if (RecoDecay::sqrtSumOfSquares(candidate.xdecayvtxv0(), candidate.ydecayvtxv0()) < V0Radius) {
        continue;
      }
      // cosPA (LFcut)
      if (candidate.cospacasc() < CascCosPA) {
        continue;
      }
      if (candidate.cospav0() < V0CosPA) {
        continue;
      }
      // cascade and v0 daughters dca cut (LF cut)
      if (candidate.dcacascdau() > dcaCascDau) {
        continue;
      }
      if (candidate.dcav0dau() > dcaV0Dau) {
        continue;
      }

      // dca omegac daughters cut
      if (candidate.dcaomegacdau() > dcaOmegacDau) {
        continue;
      }

      // cut on primary pion dcaXY
      if ((candidate.dcaxytopvprimarypi() < PrimPiMindcaxyMin) || (candidate.dcaxytopvprimarypi() > PrimPiMindcaxyMax)) {
        continue;
      }

      // pT selections
      double ptpifromcasc = sqrt((candidate.pxpifromcascatprod() * candidate.pxpifromcascatprod()) + (candidate.pypifromcascatprod() * candidate.pypifromcascatprod()));
      double ptpifromome = sqrt((candidate.pxprimarypiatprod() * candidate.pxprimarypiatprod()) + (candidate.pyprimarypiatprod() * candidate.pyprimarypiatprod()));
      if (abs(ptpifromcasc) > pTMinPiFromCasc) {
        continue;
      }
      if (abs(ptpifromome) > pTMinPiFromOme) {
        continue;
      }

      // TPC clusters selections (see O2Physics/DPG/Tasks/AOTTrack/qaEventTrack.h)
      if (trackPiFromOmeg.tpcNClsFound() < tpcClusters) { // Clusters found in TPC
        continue;
      }
      if (trackPiFromLam.tpcNClsFound() < tpcClusters) { // LF cut in cascadeanalysis.cxx
        continue;
      }
      if (trackPrFromLam.tpcNClsFound() < tpcClusters) { // LF cut in cascadeanalysis.cxx
        continue;
      }
      if (trackPiFromCasc.tpcNClsFound() < tpcClusters) { // LF cut in cascadeanalysis.cxx
        continue;
      }

      if (trackPiFromOmeg.tpcNClsCrossedRows() < tpcCrossedRows) { // Crossed rows found in TPC
        continue;
      }
      if (trackPiFromLam.tpcNClsCrossedRows() < tpcCrossedRows) {
        continue;
      }
      if (trackPrFromLam.tpcNClsCrossedRows() < tpcCrossedRows) {
        continue;
      }
      if (trackPiFromCasc.tpcNClsCrossedRows() < tpcCrossedRows) {
        continue;
      }

      if (trackPiFromOmeg.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatio) { // Crossed rows over findable clusters in TPC
        continue;
      }
      /*if(trackPiFromLam.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatio){
        continue;
      }
      if(trackPrFromLam.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatio){
        continue;
      }
      if(trackPiFromCasc.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatio){
        continue;
      }*/

      // ITS clusters selection (see O2Physics/DPG/Tasks/AOTTrack/qaEventTrack.h)
      if (trackPiFromOmeg.itsNCls() < itsClusters) { // Clusters found in ITS
        continue;
      }
      if (trackPiFromOmeg.itsNClsInnerBarrel() < itsClustersInnBarr) { // Clusters found in the inner barrel of the ITS
        continue;
      }
      if (trackPiFromCasc.itsNCls() < itsClusters) {
        continue;
      }

      // track-level PID selection
      int pidProton = -999;
      int pidPiFromLam = -999;
      int pidPiFromCasc = -999;
      int pidPiFromOme = -999;
      if (PIDTPCOnly) {
        pidProton = selectorProton.getStatusTrackPIDTPC(trackPrFromLam);
        pidPiFromLam = selectorPionFromV0.getStatusTrackPIDTPC(trackPiFromLam);
        pidPiFromCasc = selectorPionFromCasc.getStatusTrackPIDTPC(trackPiFromCasc);
        pidPiFromOme = selectorPionFromOme.getStatusTrackPIDTPC(trackPiFromOmeg);
      } else if (PIDTPCTOFCombined) {
        pidProton = selectorProton.getStatusTrackPIDAll(trackPrFromLam);
        pidPiFromLam = selectorPionFromV0.getStatusTrackPIDAll(trackPiFromLam);
        pidPiFromCasc = selectorPionFromCasc.getStatusTrackPIDAll(trackPiFromCasc);
        pidPiFromOme = selectorPionFromOme.getStatusTrackPIDAll(trackPiFromOmeg);
      }

      int statuspidLambda = -1;
      int statuspidCascade = -1;
      int statuspidOmegac = -1;

      if (pidProton == TrackSelectorPID::Status::PIDAccepted && pidPiFromLam == TrackSelectorPID::Status::PIDAccepted) {
        statuspidLambda = 1;
        hTest2->Fill(0.5);
      }

      if (pidProton == TrackSelectorPID::Status::PIDAccepted && pidPiFromLam == TrackSelectorPID::Status::PIDAccepted && pidPiFromCasc == TrackSelectorPID::Status::PIDAccepted) {
        statuspidCascade = 1;
        hTest2->Fill(1.5);
      }

      if (pidProton == TrackSelectorPID::Status::PIDAccepted && pidPiFromLam == TrackSelectorPID::Status::PIDAccepted && pidPiFromCasc == TrackSelectorPID::Status::PIDAccepted && pidPiFromOme == TrackSelectorPID::Status::PIDAccepted) {
        statuspidOmegac = 1;
        hTest2->Fill(2.5);
      }

      // invariant mass cuts
      int statusinvmassLambda = -1;
      int statusinvmassCascade = -1;
      int statusinvmassOmegac = -1;

      double invmasslambda = 0;
      if (signdecay < 0) {
        invmasslambda = candidate.invmasslambda();
      } else if (signdecay > 0) {
        invmasslambda = candidate.invmassantilambda();
      }
      double invmasscascade = candidate.invmasscascade();
      double invmassomegac = candidate.invmassomegac();

      if (abs(invmasslambda - 1.11568) < (nsigmainvmasscut * cutinvmasslambda)) {
        statusinvmassLambda = 1;
        if (statuspidLambda == 1 && statuspidCascade == 1 && statuspidOmegac == 1) {
          hTest2->Fill(3.5);
        }
      }

      if (abs(invmasscascade - 1.32171) < (nsigmainvmasscut * cutinvmasscascade)) {
        statusinvmassCascade = 1;
        if (statuspidLambda == 1 && statuspidCascade == 1 && statuspidOmegac == 1 && statusinvmassLambda == 1) {
          hTest2->Fill(4.5);
        }
      }

      // if (abs(invmassomegac-2.6952) < rangeinvmassomegac){
      if ((invmassomegac >= LowerLimitSpectrum) && (invmassomegac <= UpperLimitSpectrum)) {
        statusinvmassOmegac = 1;
        if (statuspidLambda == 1 && statuspidCascade == 1 && statuspidOmegac == 1 && statusinvmassLambda == 1 && statusinvmassCascade == 1) {
          hTest2->Fill(5.5);
        };
      }

      hfSelOmegacCandidate(statuspidLambda, statuspidCascade, statuspidOmegac, statusinvmassLambda, statusinvmassCascade, statusinvmassOmegac);

      if (statuspidLambda == -1) {
        hTest1->Fill(0.5);
      }
      if (statuspidLambda == 1) {
        hTest1->Fill(1.5);
      }
      if (statuspidCascade == -1) {
        hTest1->Fill(2.5);
      }
      if (statuspidCascade == 1) {
        hTest1->Fill(3.5);
      }
      if (statuspidOmegac == -1) {
        hTest1->Fill(4.5);
      }
      if (statuspidOmegac == 1) {
        hTest1->Fill(5.5);
      }
      if (statusinvmassLambda == -1) {
        hTest1->Fill(6.5);
      }
      if (statusinvmassLambda == 1) {
        hTest1->Fill(7.5);
      }
      if (statusinvmassCascade == -1) {
        hTest1->Fill(8.5);
      }
      if (statusinvmassCascade == 1) {
        hTest1->Fill(9.5);
      }
      if (statusinvmassOmegac == -1) {
        hTest1->Fill(10.5);
      }
      if (statusinvmassOmegac == 1) {
        hTest1->Fill(11.5);
      }

      double ptprimarypi = sqrt((candidate.pxprimarypiatprod() * candidate.pxprimarypiatprod()) + (candidate.pyprimarypiatprod() * candidate.pyprimarypiatprod()));
      if (statuspidLambda == 1 && statuspidCascade == 1 && statuspidOmegac == 1 && statusinvmassLambda == 1 && statusinvmassCascade == 1 && statusinvmassOmegac == 1) {
        // if(statuspidLambda == 1 && statuspidCascade== 1 && statuspidOmegac == 1  && statusinvmassOmegac == 1){
        hPtPrimaryPi->Fill(ptprimarypi);
        hxVertexOmegac->Fill(candidate.xdecayvtxomegac());
        hInvMassOmegac->Fill(invmassomegac);
        hCTauOmegac->Fill(candidate.ctauomegac());
        hInvMassOmegacNotFixed->Fill(candidate.massomegacnotfixed());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFOmegacCandidateSelector>(cfgc, TaskName{"hf-omegac-candidate-selector"})};
}
