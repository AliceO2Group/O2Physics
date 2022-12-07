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
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University & GSI

#include "Common/Core/TrackSelectorPID.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

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
  Configurable<double> CascRadiusMin{"CascRadiusMin", 0.5, "Min cascade radius"};
  Configurable<double> V0RadiusMin{"V0RadiusMin", 2, "Min V0 radius"};
  Configurable<double> V0CosPACut{"V0CosPACut", 0.95, "Cos PA V0 cut"};
  Configurable<double> CascCosPACut{"CascCosPACut", 0.95, "Cos PA cascade cut"};
  Configurable<double> dcaCascDauMax{"dcaCascDauMax", 5.0, "Max DCA cascade daughters"};
  Configurable<double> dcaV0DauMax{"dcaV0DauMax", 5.0, "Max DCA V0 daughters"};
  Configurable<double> dcaOmegacDauMax{"dcaOmegacDauMax", 5.0, "Max DCA omegac daughters"};

  // limit charm baryon invariant mass spectrum
  Configurable<double> LowerLimitSpectrum{"LowerLimitSpectrum", 2.4, "Lower limit invariant mass spectrum charm baryon"};
  Configurable<double> UpperLimitSpectrum{"UpperLimitSpectrum", 3.0, "Upper limit invariant mass spectrum charm baryon"};

  // kinematic selections
  Configurable<double> EtaTrackMax{"EtaTrackMax", 0.8, "Max absolute value of eta"};
  Configurable<double> ptMinPiFromCasc{"ptMinPiFromCasc", 0.15, "Min pT pi <- casc"};
  Configurable<double> ptMinPiFromOme{"ptMinPiFromOme", 0.2, "Min pT pi <- omegac"};
  Configurable<double> dcaxyPriPiMin{"dcaxyPriPiMin", 0., "Min dcaxy primary pi track to PV"};
  Configurable<double> dcaxyPriPiMax{"dcaxyPriPiMax", 10., "Max dcaxy primary pi track to PV"};

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};

  // PID options
  Configurable<bool> PidTpcOnly{"PidTpcOnly", true, "PID with TPC only"};
  Configurable<bool> PidTpcTofCombined{"PidTpcTofCombined", false, "PID with TPC & TOF"};

  // PID - TPC selections
  Configurable<double> pidTpcPtMin{"pidTpcPtMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> pidTpcPtMax{"pidTpcPtMax", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpc{"nSigmaTpc", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombined{"nSigmaTpcCombined", 5., "Nsigma cut on TPC combined with TOF"};
  // PID - TOF selections
  Configurable<double> pidTofPtMin{"pidTofPtMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> pidTofPtMax{"pidTofPtMax", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTof{"nSigmaTof", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombined{"nSigmaTofCombined", 5., "Nsigma cut on TOF combined with TPC"};

  // invariant mass cuts
  Configurable<double> SigmaInvMassLambda{"SigmaInvMassLambda", 0.0025, "Invariant mass cut for lambda (sigma)"};
  Configurable<double> SigmaInvMassCascade{"SigmaInvMassCascade", 0.0025, "Invariant mass cut for cascade (sigma)"};
  Configurable<int> nSigmaInvMassCut{"nSigmaInvMassCut", 4, "Number of sigma for invariant mass cut"};

  // detector clusters selections
  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> TpcCrossedRowsMin{"TpcCrossedRowsMin", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> TpcCrossedRowsOverFindableClustersRatioMin{"TpcCrossedRowsOverFindableClustersRatioMin", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pi <- Omegac"};
  Configurable<int> nClustersItsInnBarrMin{"nClustersItsInnBarrMin", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- Omegac"};

  using MyTrackInfo = soa::Join<aod::BigTracksPIDExtended, aod::TrackSelection>;

  OutputObj<TH1F> hPtPrimaryPi{TH1F("hPtPrimaryPi", "p_T primary #pi;p_T (GeV/#it{c});entries", 500, 0, 20)};
  OutputObj<TH1F> hxVertexOmegac{TH1F("hxVertexOmegac", "x Omegac vertex;xVtx;entries", 500, -10, 10)};
  OutputObj<TH1F> hInvMassOmegac{TH1F("hInvMassOmegac", "Omegac invariant mass;inv mass;entries", 500, 2.2, 3.1)};
  OutputObj<TH1F> hCTauOmegac{TH1F("hCTauOmegac", "Omegac ctau;ctau;entries", 500, 0., 10.)};
  OutputObj<TH1F> hInvMassOmegacNotFixed{TH1F("hInvMassOmegacNotFixed", "Omegac invariant mass (not fixed);inv mass;entries", 500, 2.2, 3.1)};

  OutputObj<TH1F> hTest1{TH1F("hTest1", "Test status steps;status;entries", 12, 0., 12.)};
  OutputObj<TH1F> hTest2{TH1F("hTest2", "Test status consecutive;status;entries", 12, 0., 12.)};

  void process(aod::HfCandOmegacBase const& candidates, MyTrackInfo const&)
  {
    TrackSelectorPID selectorPionFromOme(kPiPlus);
    selectorPionFromOme.setRangePtTPC(pidTpcPtMin, pidTpcPtMax);
    selectorPionFromOme.setRangeNSigmaTPC(-nSigmaTpc, nSigmaTpc);
    selectorPionFromOme.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombined, nSigmaTpcCombined);
    selectorPionFromOme.setRangePtTOF(pidTofPtMin, pidTofPtMax);
    selectorPionFromOme.setRangeNSigmaTOF(-nSigmaTof, nSigmaTof);
    selectorPionFromOme.setRangeNSigmaTOFCondTPC(-nSigmaTofCombined, nSigmaTofCombined);

    TrackSelectorPID selectorPionFromCasc(kPiPlus);
    selectorPionFromCasc.setRangePtTPC(pidTpcPtMin, pidTpcPtMax);
    selectorPionFromCasc.setRangeNSigmaTPC(-nSigmaTpc, nSigmaTpc);
    selectorPionFromCasc.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombined, nSigmaTpcCombined);
    selectorPionFromCasc.setRangePtTOF(pidTofPtMin, pidTofPtMax);
    selectorPionFromCasc.setRangeNSigmaTOF(-nSigmaTof, nSigmaTof);
    selectorPionFromCasc.setRangeNSigmaTOFCondTPC(-nSigmaTofCombined, nSigmaTofCombined);

    TrackSelectorPID selectorPionFromV0(kPiPlus);
    selectorPionFromV0.setRangePtTPC(pidTpcPtMin, pidTpcPtMax);
    selectorPionFromV0.setRangeNSigmaTPC(-nSigmaTpc, nSigmaTpc);
    selectorPionFromV0.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombined, nSigmaTpcCombined);
    selectorPionFromV0.setRangePtTOF(pidTofPtMin, pidTofPtMax);
    selectorPionFromV0.setRangeNSigmaTOF(-nSigmaTof, nSigmaTof);
    selectorPionFromV0.setRangeNSigmaTOFCondTPC(-nSigmaTofCombined, nSigmaTofCombined);

    TrackSelectorPID selectorProton(kProton);
    selectorProton.setRangePtTPC(pidTpcPtMin, pidTpcPtMax);
    selectorProton.setRangeNSigmaTPC(-nSigmaTpc, nSigmaTpc);
    selectorProton.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombined, nSigmaTpcCombined);
    selectorProton.setRangePtTOF(pidTofPtMin, pidTofPtMax);
    selectorProton.setRangeNSigmaTOF(-nSigmaTof, nSigmaTof);
    selectorProton.setRangeNSigmaTOFCondTPC(-nSigmaTofCombined, nSigmaTofCombined);

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
      if (abs(etav0posdau) > EtaTrackMax) {
        continue;
      }
      if (abs(etav0negdau) > EtaTrackMax) {
        continue;
      }
      if (abs(etapifromcasc) > EtaTrackMax) {
        continue;
      }
      if (abs(etapifromome) > EtaTrackMax) {
        continue;
      }

      // minimum radius cut (LFcut)
      if (RecoDecay::sqrtSumOfSquares(candidate.xdecayvtxcascade(), candidate.ydecayvtxcascade()) < CascRadiusMin) {
        continue;
      }
      if (RecoDecay::sqrtSumOfSquares(candidate.xdecayvtxv0(), candidate.ydecayvtxv0()) < V0RadiusMin) {
        continue;
      }
      // cosPA (LFcut)
      if (candidate.cospacasc() < CascCosPACut) {
        continue;
      }
      if (candidate.cospav0() < V0CosPACut) {
        continue;
      }
      // cascade and v0 daughters dca cut (LF cut)
      if (candidate.dcacascdau() > dcaCascDauMax) {
        continue;
      }
      if (candidate.dcav0dau() > dcaV0DauMax) {
        continue;
      }

      // dca omegac daughters cut
      if (candidate.dcaomegacdau() > dcaOmegacDauMax) {
        continue;
      }

      // cut on primary pion dcaXY
      if ((candidate.dcaxytopvprimarypi() < dcaxyPriPiMin) || (candidate.dcaxytopvprimarypi() > dcaxyPriPiMax)) {
        continue;
      }

      // pT selections
      double ptpifromcasc = sqrt((candidate.pxpifromcascatprod() * candidate.pxpifromcascatprod()) + (candidate.pypifromcascatprod() * candidate.pypifromcascatprod()));
      double ptpifromome = sqrt((candidate.pxprimarypiatprod() * candidate.pxprimarypiatprod()) + (candidate.pyprimarypiatprod() * candidate.pyprimarypiatprod()));
      if (abs(ptpifromcasc) > ptMinPiFromCasc) {
        continue;
      }
      if (abs(ptpifromome) > ptMinPiFromOme) {
        continue;
      }

      //  TPC clusters selections
      if (trackPiFromOmeg.tpcNClsFound() < nClustersTpcMin) {
        continue;
      }
      if (trackPiFromLam.tpcNClsFound() < nClustersTpcMin) {
        continue;
      }
      if (trackPrFromLam.tpcNClsFound() < nClustersTpcMin) {
        continue;
      }
      if (trackPiFromCasc.tpcNClsFound() < nClustersTpcMin) {
        continue;
      }

      if (trackPiFromOmeg.tpcNClsCrossedRows() < TpcCrossedRowsMin) {
        continue;
      }
      if (trackPiFromLam.tpcNClsCrossedRows() < TpcCrossedRowsMin) {
        continue;
      }
      if (trackPrFromLam.tpcNClsCrossedRows() < TpcCrossedRowsMin) {
        continue;
      }
      if (trackPiFromCasc.tpcNClsCrossedRows() < TpcCrossedRowsMin) {
        continue;
      }

      if (trackPiFromOmeg.tpcCrossedRowsOverFindableCls() < TpcCrossedRowsOverFindableClustersRatioMin) {
        continue;
      }

      //  ITS clusters selection
      if (trackPiFromOmeg.itsNCls() < nClustersItsMin) {
        continue;
      }
      if (trackPiFromOmeg.itsNClsInnerBarrel() < nClustersItsInnBarrMin) {
        continue;
      }
      if (trackPiFromCasc.itsNCls() < nClustersItsMin) {
        continue;
      }

      // track-level PID selection
      int pidProton = -999;
      int pidPiFromLam = -999;
      int pidPiFromCasc = -999;
      int pidPiFromOme = -999;
      if (PidTpcOnly) {
        pidProton = selectorProton.getStatusTrackPIDTPC(trackPrFromLam);
        pidPiFromLam = selectorPionFromV0.getStatusTrackPIDTPC(trackPiFromLam);
        pidPiFromCasc = selectorPionFromCasc.getStatusTrackPIDTPC(trackPiFromCasc);
        pidPiFromOme = selectorPionFromOme.getStatusTrackPIDTPC(trackPiFromOmeg);
      } else if (PidTpcTofCombined) {
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

      if (abs(invmasslambda - 1.11568) < (nSigmaInvMassCut * SigmaInvMassLambda)) {
        statusinvmassLambda = 1;
        if (statuspidLambda == 1 && statuspidCascade == 1 && statuspidOmegac == 1) {
          hTest2->Fill(3.5);
        }
      }

      if (abs(invmasscascade - 1.32171) < (nSigmaInvMassCut * SigmaInvMassCascade)) {
        statusinvmassCascade = 1;
        if (statuspidLambda == 1 && statuspidCascade == 1 && statuspidOmegac == 1 && statusinvmassLambda == 1) {
          hTest2->Fill(4.5);
        }
      }

      if ((invmassomegac >= LowerLimitSpectrum) && (invmassomegac <= UpperLimitSpectrum)) {
        statusinvmassOmegac = 1;
        if (statuspidLambda == 1 && statuspidCascade == 1 && statuspidOmegac == 1 && statusinvmassLambda == 1 && statusinvmassCascade == 1) {
          hTest2->Fill(5.5);
        }
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
