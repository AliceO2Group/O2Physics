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

/// \file candidateSelectorOmegac.cxx
/// \brief Omegac → Xi Pi selection task
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University & GSI

#include "Common/Core/TrackSelectorPID.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/RecoDecay.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis::pdg;
using namespace o2::aod::hf_cand_toxipi;
using namespace o2::aod::hf_sel_toxipi;

/// Struct for applying Omegac selection cuts
struct HfCandidateSelectorToXiPi {
  Produces<aod::HfSelToXiPi> hfSelToXiPi;

  // LF analysis selections
  // zPV -> can be already set in HFeventselection -> 10 cm
  // sel8 -> can be already set in HFeventselection -> true
  Configurable<double> radiusCascMin{"radiusCascMin", 0.5, "Min cascade radius"};
  Configurable<double> radiusV0Min{"radiusV0Min", 2, "Min V0 radius"};
  Configurable<double> cosPAV0Min{"cosPAV0Min", 0.95, "Min valueCosPA V0"};
  Configurable<double> cosPACascMin{"cosPACascMin", 0.95, "Min value CosPA cascade"};
  Configurable<double> dcaCascDauMax{"dcaCascDauMax", 5.0, "Max DCA cascade daughters"};
  Configurable<double> dcaV0DauMax{"dcaV0DauMax", 5.0, "Max DCA V0 daughters"};
  Configurable<double> dcaOmegacDauMax{"dcaOmegacDauMax", 5.0, "Max DCA omegac daughters"};

  // limit charm baryon invariant mass spectrum
  Configurable<double> invMassOmegacMin{"invMassOmegacMin", 2.0, "Lower limit invariant mass spectrum charm baryon"}; // 2.4 Omegac0 only
  Configurable<double> invMassOmegacMax{"invMassOmegacMax", 3.1, "Upper limit invariant mass spectrum charm baryon"};

  // kinematic selections
  Configurable<double> etaTrackMax{"etaTrackMax", 0.8, "Max absolute value of eta"};
  Configurable<double> ptPiFromCascMin{"ptPiFromCascMin", 0.15, "Min pT pi <- casc"};
  Configurable<double> ptPiFromOmeMin{"ptPiFromOmeMin", 0.2, "Min pT pi <- omegac"};
  Configurable<double> dcaXYPriPiMin{"dcaXYPriPiMin", 0., "Min dcaxy primary pi track to PV"};
  Configurable<double> dcaXYPriPiMax{"dcaXYPriPiMax", 10., "Max dcaxy primary pi track to PV"};

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};

  // PID options
  Configurable<bool> usePidTpcOnly{"usePidTpcOnly", true, "Perform PID using only TPC"};
  Configurable<bool> usePidTpcTofCombined{"usePidTpcTofCombined", false, "Perform PID using TPC & TOF"};

  // PID - TPC selections
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};

  // PID - TOF selections
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};

  // invariant mass cuts
  Configurable<double> sigmaInvMassLambda{"sigmaInvMassLambda", 0.0025, "Invariant mass cut for lambda (sigma)"};
  Configurable<double> sigmaInvMassCascade{"sigmaInvMassCascade", 0.0025, "Invariant mass cut for cascade (sigma)"};
  Configurable<int> nSigmaInvMassCut{"nSigmaInvMassCut", 4, "Number of sigma for invariant mass cut"};

  // detector clusters selections
  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> nTpcCrossedRowsMin{"nTpcCrossedRowsMin", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> tpcCrossedRowsOverFindableClustersRatioMin{"tpcCrossedRowsOverFindableClustersRatioMin", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pi <- Omegac"};
  Configurable<int> nClustersItsInnBarrMin{"nClustersItsInnBarrMin", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- Omegac"};

  using MyTrackInfo = soa::Join<aod::BigTracksPIDExtended, aod::TrackSelection>;

  OutputObj<TH1F> hxVertexOmegac{TH1F("hxVertexOmegac", "x Omegac vertex;xVtx;entries", 500, -10, 10)};
  OutputObj<TH1F> hInvMassOmegac{TH1F("hInvMassOmegac", "Omegac invariant mass;inv mass;entries", 500, 2.2, 3.1)};
  OutputObj<TH1F> hCTauOmegac{TH1F("hCTauOmegac", "Omegac ctau;ctau;entries", 500, 0., 10.)};
  OutputObj<TH1F> hCTauXic{TH1F("hCTauXic", "Xic ctau;ctau;entries", 500, 0., 10.)};

  // temporary histo for debugging (to be removed after test on hyperloop)
  OutputObj<TH1F> hTest1{TH1F("hTest1", "Test status steps;status;entries", 12, 0., 12.)};
  OutputObj<TH1F> hTest2{TH1F("hTest2", "Test status consecutive;status;entries", 12, 0., 12.)};

  void process(aod::HfCandToXiPi const& candidates, MyTrackInfo const&)
  {
    TrackSelectorPID selectorPionFromOme(kPiPlus);
    selectorPionFromOme.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
    selectorPionFromOme.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPionFromOme.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPionFromOme.setRangePtTOF(ptPidTofMin, ptPidTofMax);
    selectorPionFromOme.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);
    selectorPionFromOme.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);

    TrackSelectorPID selectorPionFromCasc(kPiPlus);
    selectorPionFromCasc.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
    selectorPionFromCasc.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPionFromCasc.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPionFromCasc.setRangePtTOF(ptPidTofMin, ptPidTofMax);
    selectorPionFromCasc.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);
    selectorPionFromCasc.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);

    TrackSelectorPID selectorPionFromV0(kPiPlus);
    selectorPionFromV0.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
    selectorPionFromV0.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPionFromV0.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPionFromV0.setRangePtTOF(ptPidTofMin, ptPidTofMax);
    selectorPionFromV0.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);
    selectorPionFromV0.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);

    TrackSelectorPID selectorProton(kProton);
    selectorProton.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
    selectorProton.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
    selectorProton.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorProton.setRangePtTOF(ptPidTofMin, ptPidTofMax);
    selectorProton.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);
    selectorProton.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);

    double massLambdaFromPDG = RecoDecay::getMassPDG(kLambda0);
    double massXiFromPDG = RecoDecay::getMassPDG(kXiMinus);

    // looping over omegac candidates
    for (auto const& candidate : candidates) {

      bool resultSelections = true; // True if the candidate passes all the selections, False otherwise

      auto trackV0PosDau = candidate.posTrack_as<MyTrackInfo>();    // positive V0 daughter
      auto trackV0NegDau = candidate.negTrack_as<MyTrackInfo>();    // negative V0 daughter
      auto trackPiFromCasc = candidate.bachelor_as<MyTrackInfo>();  // pion <- cascade
      auto trackPiFromOmeg = candidate.primaryPi_as<MyTrackInfo>(); // pion <- omegac

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      int8_t signDecay = candidate.signDecay(); // sign of pi <- cascade

      if (signDecay > 0) {
        trackPiFromLam = trackV0PosDau;
        trackPrFromLam = trackV0NegDau;
      } else if (signDecay == 0) {
        resultSelections = false;
      }

      // eta selection
      double etaV0PosDau = candidate.etaV0PosDau();
      double etaV0NegDau = candidate.etaV0NegDau();
      double etaPiFromCasc = candidate.etaPiFromCasc();
      double etaPiFromOme = candidate.etaPiFromOme();
      if (std::abs(etaV0PosDau) > etaTrackMax) {
        resultSelections = false;
      }
      if (std::abs(etaV0NegDau) > etaTrackMax) {
        resultSelections = false;
      }
      if (std::abs(etaPiFromCasc) > etaTrackMax) {
        resultSelections = false;
      }
      if (std::abs(etaPiFromOme) > etaTrackMax) {
        resultSelections = false;
      }

      // minimum radius cut (LFcut)
      if (RecoDecay::sqrtSumOfSquares(candidate.xDecayVtxCascade(), candidate.yDecayVtxCascade()) < radiusCascMin) {
        resultSelections = false;
      }
      if (RecoDecay::sqrtSumOfSquares(candidate.xDecayVtxV0(), candidate.yDecayVtxV0()) < radiusV0Min) {
        resultSelections = false;
      }
      // cosPA (LFcut)
      if (candidate.cosPACasc() < cosPACascMin) {
        resultSelections = false;
      }
      if (candidate.cosPAV0() < cosPAV0Min) {
        resultSelections = false;
      }
      // cascade and v0 daughters dca cut (LF cut)
      if (candidate.dcaCascDau() > dcaCascDauMax) {
        resultSelections = false;
      }
      if (candidate.dcaV0Dau() > dcaV0DauMax) {
        resultSelections = false;
      }

      // dca omegac daughters cut
      if (candidate.dcaOmegacDau() > dcaOmegacDauMax) {
        resultSelections = false;
      }

      // cut on primary pion dcaXY
      if ((candidate.dcaXYToPVPrimaryPi() < dcaXYPriPiMin) || (candidate.dcaXYToPVPrimaryPi() > dcaXYPriPiMax)) {
        resultSelections = false;
      }

      // pT selections
      double ptPiFromCasc = RecoDecay::sqrtSumOfSquares(candidate.pxPiFromCasc(), candidate.pyPiFromCasc());
      double ptPiFromOme = RecoDecay::sqrtSumOfSquares(candidate.pxPrimaryPi(), candidate.pyPrimaryPi());
      if (std::abs(ptPiFromCasc) < ptPiFromCascMin) {
        resultSelections = false;
      }
      if (std::abs(ptPiFromOme) < ptPiFromOmeMin) {
        resultSelections = false;
      }

      //  TPC clusters selections
      if (trackPiFromOmeg.tpcNClsFound() < nClustersTpcMin) {
        resultSelections = false;
      }
      if (trackPiFromLam.tpcNClsFound() < nClustersTpcMin) {
        resultSelections = false;
      }
      if (trackPrFromLam.tpcNClsFound() < nClustersTpcMin) {
        resultSelections = false;
      }
      if (trackPiFromCasc.tpcNClsFound() < nClustersTpcMin) {
        resultSelections = false;
      }

      if (trackPiFromOmeg.tpcNClsCrossedRows() < nTpcCrossedRowsMin) {
        resultSelections = false;
      }
      if (trackPiFromLam.tpcNClsCrossedRows() < nTpcCrossedRowsMin) {
        resultSelections = false;
      }
      if (trackPrFromLam.tpcNClsCrossedRows() < nTpcCrossedRowsMin) {
        resultSelections = false;
      }
      if (trackPiFromCasc.tpcNClsCrossedRows() < nTpcCrossedRowsMin) {
        resultSelections = false;
      }

      if (trackPiFromOmeg.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatioMin) {
        resultSelections = false;
      }

      //  ITS clusters selection
      if (trackPiFromOmeg.itsNCls() < nClustersItsMin) {
        resultSelections = false;
      }
      if (trackPiFromOmeg.itsNClsInnerBarrel() < nClustersItsInnBarrMin) {
        resultSelections = false;
      }
      if (trackPiFromCasc.itsNCls() < nClustersItsMin) {
        resultSelections = false;
      }

      // track-level PID selection
      int pidProton = -999;
      int pidPiFromLam = -999;
      int pidPiFromCasc = -999;
      int pidPiFromOme = -999;
      if (usePidTpcOnly == usePidTpcTofCombined) {
        LOGF(fatal, "Check the PID configurables, usePidTpcOnly and usePidTpcTofCombined can't have the same value");
      }
      if (usePidTpcOnly) {
        pidProton = selectorProton.getStatusTrackPIDTPC(trackPrFromLam);
        pidPiFromLam = selectorPionFromV0.getStatusTrackPIDTPC(trackPiFromLam);
        pidPiFromCasc = selectorPionFromCasc.getStatusTrackPIDTPC(trackPiFromCasc);
        pidPiFromOme = selectorPionFromOme.getStatusTrackPIDTPC(trackPiFromOmeg);
      } else if (usePidTpcTofCombined) {
        pidProton = selectorProton.getStatusTrackPIDAll(trackPrFromLam);
        pidPiFromLam = selectorPionFromV0.getStatusTrackPIDAll(trackPiFromLam);
        pidPiFromCasc = selectorPionFromCasc.getStatusTrackPIDAll(trackPiFromCasc);
        pidPiFromOme = selectorPionFromOme.getStatusTrackPIDAll(trackPiFromOmeg);
      }

      bool statusPidLambda = false;
      bool statusPidCascade = false;
      bool statusPidOmegac = false;

      if (pidProton == TrackSelectorPID::Status::PIDAccepted && pidPiFromLam == TrackSelectorPID::Status::PIDAccepted) {
        statusPidLambda = true;
        if (resultSelections) {
          hTest2->Fill(0.5);
        }
      }

      if (pidProton == TrackSelectorPID::Status::PIDAccepted && pidPiFromLam == TrackSelectorPID::Status::PIDAccepted && pidPiFromCasc == TrackSelectorPID::Status::PIDAccepted) {
        statusPidCascade = true;
        if (resultSelections) {
          hTest2->Fill(1.5);
        }
      }

      if (pidProton == TrackSelectorPID::Status::PIDAccepted && pidPiFromLam == TrackSelectorPID::Status::PIDAccepted && pidPiFromCasc == TrackSelectorPID::Status::PIDAccepted && pidPiFromOme == TrackSelectorPID::Status::PIDAccepted) {
        statusPidOmegac = true;
        if (resultSelections) {
          hTest2->Fill(2.5);
        }
      }

      // invariant mass cuts
      bool statusInvMassLambda = false;
      bool statusInvMassCascade = false;
      bool statusInvMassOmegac = false;

      double invMassLambda = 0;
      if (signDecay < 0) {
        invMassLambda = candidate.invMassLambda();
      } else if (signDecay > 0) {
        invMassLambda = candidate.invMassAntiLambda();
      }
      double invMassCascade = candidate.invMassCascade();
      double invMassOmegac = candidate.invMassOmegac();

      if (std::abs(invMassLambda - massLambdaFromPDG) < (nSigmaInvMassCut * sigmaInvMassLambda)) {
        statusInvMassLambda = true;
        if (statusPidLambda && statusPidCascade && statusPidOmegac && resultSelections) {
          hTest2->Fill(3.5);
        }
      }

      if (std::abs(invMassCascade - massXiFromPDG) < (nSigmaInvMassCut * sigmaInvMassCascade)) {
        statusInvMassCascade = true;
        if (statusPidLambda && statusPidCascade && statusPidOmegac && statusInvMassLambda && resultSelections) {
          hTest2->Fill(4.5);
        }
      }

      if ((invMassOmegac >= invMassOmegacMin) && (invMassOmegac <= invMassOmegacMax)) {
        statusInvMassOmegac = true;
        if (statusPidLambda && statusPidCascade && statusPidOmegac && statusInvMassLambda && statusInvMassCascade && resultSelections) {
          hTest2->Fill(5.5);
        }
      }

      hfSelToXiPi(statusPidLambda, statusPidCascade, statusPidOmegac, statusInvMassLambda, statusInvMassCascade, statusInvMassOmegac, resultSelections);

      if (resultSelections) {
        if (!statusPidLambda) {
          hTest1->Fill(0.5);
        }
        if (statusPidLambda) {
          hTest1->Fill(1.5);
        }
        if (!statusPidCascade) {
          hTest1->Fill(2.5);
        }
        if (statusPidCascade) {
          hTest1->Fill(3.5);
        }
        if (!statusPidOmegac) {
          hTest1->Fill(4.5);
        }
        if (statusPidOmegac) {
          hTest1->Fill(5.5);
        }
        if (!statusInvMassLambda) {
          hTest1->Fill(6.5);
        }
        if (statusInvMassLambda) {
          hTest1->Fill(7.5);
        }
        if (!statusInvMassCascade) {
          hTest1->Fill(8.5);
        }
        if (statusInvMassCascade) {
          hTest1->Fill(9.5);
        }
        if (!statusInvMassOmegac) {
          hTest1->Fill(10.5);
        }
        if (statusInvMassOmegac) {
          hTest1->Fill(11.5);
        }
      }

      if (statusPidLambda && statusPidCascade && statusPidOmegac && statusInvMassLambda && statusInvMassCascade && statusInvMassOmegac && resultSelections) {
        hxVertexOmegac->Fill(candidate.xDecayVtxOmegac());
        hInvMassOmegac->Fill(invMassOmegac);
        hCTauOmegac->Fill(candidate.ctauOmegac());
        hCTauXic->Fill(candidate.ctauXic());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorToXiPi>(cfgc)};
}
