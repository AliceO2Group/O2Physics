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

/// \file taskDirectedFlowCharmHadrons.cxx
/// \brief Analysis task for charm hadron directed flow
///
/// \author Prottay Das, prottay.das@cern.ch
/// \author Biao Zhang, biao.zhanng@cern.ch

#include <string>
#include <vector>
#include <TMath.h>
#include <cstdlib>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/EventPlaneHelper.h"
#include "PWGLF/DataModel/SPCalibrationTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;

enum DecayChannel { DplusToPiKPi = 0,
                    D0ToPiK,
                    D0ToKPi,
                    DstarToD0Pi };

struct HfTaskDirectedFlowCharmHadrons {
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<bool> selectionFlagDstar{"selectionFlagDstar", false, "Selection Flag for Dstar"};
  Configurable<int> selectionFlag{"selectionFlag", 1, "Selection Flag for hadron (e.g. 1 for skimming, 3 for topo. and kine., 7 for PID)"};
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality accepted in SP computation"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality accepted in SP computation"};
  Configurable<bool> storeMl{"storeMl", false, "Flag to store ML scores"};
  Configurable<bool> direct{"direct", false, "Flag to calculate direct v1 odd and even"};
  Configurable<bool> userap{"userap", false, "Flag to fill rapidity vs eta "};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indices of BDT scores to be stored. Two indexes max."};

  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {100, 1.78, 2.05}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0}, ""};
  ConfigurableAxis thnConfigAxisEta{"thnConfigAxisEta", {VARIABLE_WIDTH, -0.8, -0.4, 0, 0.4, 0.8}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0}, ""};
  ConfigurableAxis thnConfigAxisScalarProd{"thnConfigAxisScalarProd", {8000, -2.0, 2.0}, ""};
  ConfigurableAxis thnConfigAxisSign{"thnConfigAxisSign", {2, -2.0, 2.0}, ""};
  ConfigurableAxis thnConfigAxisMlOne{"thnConfigAxisMlOne", {1000, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisMlTwo{"thnConfigAxisMlTwo", {1000, 0., 1.}, ""};

  using CandDplusDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandD0DataWMl = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using CandD0Data = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  using CandDstarDataWMl = soa::Filtered<soa::Join<aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi>>;
  using CandDstarData = soa::Filtered<soa::Join<aod::HfCandDstars, aod::HfSelDstarToD0Pi>>;
  using CollsWithQvecs = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::SPCalibrationTables>;
  using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Filter filterSelectDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstar;

  Partition<CandD0Data> selectedD0ToPiK = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0Data> selectedD0ToKPi = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Partition<CandD0DataWMl> selectedD0ToPiKWMl = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0DataWMl> selectedD0ToKPiWMl = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;

  SliceCache cache;
  HfHelper hfHelper;
  EventPlaneHelper epHelper;
  HfEventSelection hfEvSel; // event selection and monitoring
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {

    /// check process functions
    std::array<int, 6> processes = {doprocessDplusStd, doprocessDplusMl, doprocessD0Std, doprocessD0Ml, doprocessDstarStd, doprocessDstarMl};
    const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    } else if (nProcesses == 0) {
      LOGP(fatal, "No process function enabled");
    }

    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "Inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisEta{thnConfigAxisEta, "#it{#eta}"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality"};
    const AxisSpec thnAxisScalarProd{thnConfigAxisScalarProd, "SP"};
    const AxisSpec thnAxisSign{thnConfigAxisSign, "Sign"};
    const AxisSpec thnAxisMlOne{thnConfigAxisMlOne, "Bkg score"};
    const AxisSpec thnAxisMlTwo{thnConfigAxisMlTwo, "FD score"};

    std::vector<AxisSpec> axes = {thnAxisInvMass, thnAxisCent, thnAxisPt, thnAxisEta, thnAxisScalarProd, thnAxisSign};
    if (storeMl) {
      axes.insert(axes.end(), {thnAxisMlOne, thnAxisMlTwo});
    }

    if (direct) {
      registry.add("hpQxytpvscent", "hpQxytpvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpuxyQxypvscentpteta", "hpuxyQxypvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpuxyQxytvscentpteta", "hpuxyQxytvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpoddvscentpteta", "hpoddvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpevenvscentpteta", "hpevenvscentpteta", HistType::kTHnSparseF, axes, true);

      registry.add("hpQxpvscent", "hpQxpvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQypvscent", "hpQypvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQxtvscent", "hpQxtvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQytvscent", "hpQytvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpuxvscentpteta", "hpuxvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpuyvscentpteta", "hpuyvscentpteta", HistType::kTHnSparseF, axes, true);
    } else {
      registry.add("hpQxtQxpvscent", "hpQxtQxpvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQytQypvscent", "hpQytQypvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQxtQypvscent", "hpQxtQypvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQxpQytvscent", "hpQxpQytvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpuxQxpvscentpteta", "hpuxQxpvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpuyQypvscentpteta", "hpuyQypvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpuxQxtvscentpteta", "hpuxQxtvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpuyQytvscentpteta", "hpuyQytvscentpteta", HistType::kTHnSparseF, axes, true);

      registry.add("hpQxpvscent", "hpQxpvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQypvscent", "hpQypvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQxtvscent", "hpQxtvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQytvscent", "hpQytvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpuxvscentpteta", "hpuxvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpuyvscentpteta", "hpuyvscentpteta", HistType::kTHnSparseF, axes, true);
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }; // end init

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  double getCentrality(CollsWithQvecs::iterator const& collision)
  {
    double cent = -999.;
    switch (centEstimator) {
      case CentralityEstimator::FV0A:
        cent = collision.centFV0A();
        break;
      case CentralityEstimator::FT0M:
        cent = collision.centFT0M();
        break;
      case CentralityEstimator::FT0A:
        cent = collision.centFT0A();
        break;
      case CentralityEstimator::FT0C:
        cent = collision.centFT0C();
        break;
      default:
        LOG(warning) << "Centrality estimator not valid. Possible values are V0A, T0M, T0A, T0C. Fallback to V0A";
        cent = collision.centFV0A();
        break;
    }
    return cent;
  }

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param candidates are the selected candidates
  template <DecayChannel channel, typename T1, typename Trk>
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision,
                       T1 const& candidates,
                       Trk const& /*tracks*/)
  {
    double cent = getCentrality(collision);
    if (cent < centralityMin || cent > centralityMax) {
      return;
    }

    if (!collision.triggereventsp()) { // for selecting only callibrated events
      return;
    }

    auto qxZDCA = collision.qxZDCA();
    auto qyZDCA = collision.qyZDCA();
    auto qxZDCC = collision.qxZDCC(); // extracting q vectors of ZDC
    auto qyZDCC = collision.qyZDCC();

    auto QxtQxp = qxZDCC * qxZDCA;
    auto QytQyp = qyZDCC * qyZDCA;
    auto Qxytp = QxtQxp + QytQyp;
    auto QxpQyt = qxZDCA * qyZDCC;
    auto QxtQyp = qxZDCC * qyZDCA;

    // correlations in the denominators for SP calculation
    if (direct) {
      registry.fill(HIST("hpQxytpvscent"), cent, Qxytp);
      registry.fill(HIST("hpQxpvscent"), cent, qxZDCA);
      registry.fill(HIST("hpQxtvscent"), cent, qxZDCC);
      registry.fill(HIST("hpQypvscent"), cent, qyZDCA);
      registry.fill(HIST("hpQytvscent"), cent, qyZDCC);
    } else {
      registry.fill(HIST("hpQxtQxpvscent"), cent, QxtQxp);
      registry.fill(HIST("hpQytQypvscent"), cent, QytQyp);
      registry.fill(HIST("hpQxpQytvscent"), cent, QxpQyt);
      registry.fill(HIST("hpQxtQypvscent"), cent, QxtQyp);
      registry.fill(HIST("hpQxpvscent"), cent, qxZDCA);
      registry.fill(HIST("hpQxtvscent"), cent, qxZDCC);
      registry.fill(HIST("hpQypvscent"), cent, qyZDCA);
      registry.fill(HIST("hpQytvscent"), cent, qyZDCC);
    }

    for (const auto& candidate : candidates) {
      double massCand = 0.;
      double rapCand = 0.;
      double sign = 0.; // electric charge of the first daughter track to differentiate particle and antiparticle
      double signDstarCand = 0.0;
      std::vector<double> outputMl = {-999., -999.};
      if constexpr (std::is_same_v<T1, CandDplusData> || std::is_same_v<T1, CandDplusDataWMl>) {
        massCand = hfHelper.invMassDplusToPiKPi(candidate);
        rapCand = hfHelper.yDplus(candidate);
        auto trackprong0 = candidate.template prong0_as<Trk>();
        sign = trackprong0.sign();
        if constexpr (std::is_same_v<T1, CandDplusDataWMl>) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
            outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
        }
      } else if constexpr (std::is_same_v<T1, CandD0Data> || std::is_same_v<T1, CandD0DataWMl>) {
        switch (channel) {
          case DecayChannel::D0ToPiK:
            massCand = hfHelper.invMassD0ToPiK(candidate);
            rapCand = hfHelper.yD0(candidate);
            sign = candidate.isSelD0bar() ? 3 : 1; // 3: reflected D0bar, 1: pure D0 excluding reflected D0bar
            if constexpr (std::is_same_v<T1, CandD0DataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
                outputMl[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
            }
            break;
          case DecayChannel::D0ToKPi:
            massCand = hfHelper.invMassD0barToKPi(candidate);
            rapCand = hfHelper.yD0(candidate);
            sign = candidate.isSelD0() ? 3 : 2; // 3: reflected D0, 2: pure D0bar excluding reflected D0
            if constexpr (std::is_same_v<T1, CandD0DataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
                outputMl[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
            }
            break;
          default:
            break;
        }
      } else if constexpr (std::is_same_v<T1, CandDstarData> || std::is_same_v<T1, CandDstarDataWMl>) {
        signDstarCand = candidate.signSoftPi();
        if (candidate.signSoftPi() > 0) {
          massCand = std::abs(candidate.invMassDstar() - candidate.invMassD0());
          rapCand = candidate.y(candidate.invMassDstar());
        } else if (candidate.signSoftPi() < 0) {
          massCand = std::abs(candidate.invMassAntiDstar() - candidate.invMassD0Bar());
          rapCand = candidate.y(candidate.invMassAntiDstar());
        }
        if constexpr (std::is_same_v<T1, CandDstarDataWMl>) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
            outputMl[iclass] = candidate.mlProbDstarToD0Pi()[classMl->at(iclass)];
        }
      }

      double ptCand = candidate.pt();
      double etaCand = candidate.eta();
      double phiCand = candidate.phi();
      double cosNPhi = std::cos(phiCand);
      double sinNPhi = std::sin(phiCand);

      if (userap)
        etaCand = rapCand;

      if (selectionFlagDstar)
        sign = signDstarCand;

      auto ux = cosNPhi; // real part of candidate q vector
      auto uy = sinNPhi; // imaginary part of candidate q vector
      auto uxQxp = ux * qxZDCA;
      auto uyQyp = uy * qyZDCA; // correlations of particle and ZDC q vectors
      auto uxyQxyp = uxQxp + uyQyp;
      auto uxQxt = ux * qxZDCC;
      auto uyQyt = uy * qyZDCC;
      auto uxyQxyt = uxQxt + uyQyt;
      auto oddv1 = ux * (qxZDCA - qxZDCC) + uy * (qyZDCA - qyZDCC);
      auto evenv1 = ux * (qxZDCA + qxZDCC) + uy * (qyZDCA + qyZDCC);

      if (storeMl) {
        if (direct) {
          registry.fill(HIST("hpuxyQxypvscentpteta"), massCand, cent, ptCand, etaCand, uxyQxyp, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpuxyQxytvscentpteta"), massCand, cent, ptCand, etaCand, uxyQxyt, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpoddvscentpteta"), massCand, cent, ptCand, etaCand, oddv1, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpevenvscentpteta"), massCand, cent, ptCand, etaCand, evenv1, sign, outputMl[0], outputMl[1]);

          registry.fill(HIST("hpuxvscentpteta"), massCand, cent, ptCand, etaCand, ux, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpuyvscentpteta"), massCand, cent, ptCand, etaCand, uy, sign, outputMl[0], outputMl[1]);

        } else {
          registry.fill(HIST("hpuxQxpvscentpteta"), massCand, cent, ptCand, etaCand, uxQxp, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpuyQypvscentpteta"), massCand, cent, ptCand, etaCand, uyQyp, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpuxQxtvscentpteta"), massCand, cent, ptCand, etaCand, uxQxt, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpuyQytvscentpteta"), massCand, cent, ptCand, etaCand, uyQyt, sign, outputMl[0], outputMl[1]);

          registry.fill(HIST("hpuxvscentpteta"), massCand, cent, ptCand, etaCand, ux, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpuyvscentpteta"), massCand, cent, ptCand, etaCand, uy, sign, outputMl[0], outputMl[1]);
        }
      } else {
        if (direct) {
          registry.fill(HIST("hpuxyQxypvscentpteta"), massCand, cent, ptCand, etaCand, uxyQxyp, sign);
          registry.fill(HIST("hpuxyQxytvscentpteta"), massCand, cent, ptCand, etaCand, uxyQxyt, sign);
          registry.fill(HIST("hpoddvscentpteta"), massCand, cent, ptCand, etaCand, oddv1, sign);
          registry.fill(HIST("hpevenvscentpteta"), massCand, cent, ptCand, etaCand, evenv1, sign);

          registry.fill(HIST("hpuxvscentpteta"), massCand, cent, ptCand, etaCand, ux, sign);
          registry.fill(HIST("hpuyvscentpteta"), massCand, cent, ptCand, etaCand, uy, sign);
        } else {
          registry.fill(HIST("hpuxQxpvscentpteta"), massCand, cent, ptCand, etaCand, uxQxp, sign);
          registry.fill(HIST("hpuyQypvscentpteta"), massCand, cent, ptCand, etaCand, uyQyp, sign);
          registry.fill(HIST("hpuxQxtvscentpteta"), massCand, cent, ptCand, etaCand, uxQxt, sign);
          registry.fill(HIST("hpuyQytvscentpteta"), massCand, cent, ptCand, etaCand, uyQyt, sign);

          registry.fill(HIST("hpuxvscentpteta"), massCand, cent, ptCand, etaCand, ux, sign);
          registry.fill(HIST("hpuyvscentpteta"), massCand, cent, ptCand, etaCand, uy, sign);
        }
      }
    }
  }
  // D0 with ML
  void processD0Ml(CollsWithQvecs::iterator const& collision,
                   CandD0DataWMl const& /*candidatesD0*/,
                   TracksWithExtra const& tracks)
  {
    auto candsD0ToPiKWMl = selectedD0ToPiKWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsD0ToKPiWMl = selectedD0ToKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::D0ToPiK>(collision, candsD0ToPiKWMl, tracks);
    runFlowAnalysis<DecayChannel::D0ToKPi>(collision, candsD0ToKPiWMl, tracks);
  }
  PROCESS_SWITCH(HfTaskDirectedFlowCharmHadrons, processD0Ml, "Process D0 candidates with ML", false);

  // D0 with rectangular cuts
  void processD0Std(CollsWithQvecs::iterator const& collision,
                    CandD0Data const& /*candidatesD0*/,
                    TracksWithExtra const& tracks)
  {
    auto candsD0ToPiK = selectedD0ToPiK->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsD0ToKPi = selectedD0ToKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::D0ToPiK>(collision, candsD0ToPiK, tracks);
    runFlowAnalysis<DecayChannel::D0ToKPi>(collision, candsD0ToKPi, tracks);
  }
  PROCESS_SWITCH(HfTaskDirectedFlowCharmHadrons, processD0Std, "Process D0 candidates with rectangular cuts", false);

  // Dplus with ML
  void processDplusMl(CollsWithQvecs::iterator const& collision,
                      CandDplusDataWMl const& candidatesDplus,
                      TracksWithExtra const& tracks)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus, tracks);
  }
  PROCESS_SWITCH(HfTaskDirectedFlowCharmHadrons, processDplusMl, "Process Dplus candidates with ML", false);

  // Dplus with rectangular cuts
  void processDplusStd(CollsWithQvecs::iterator const& collision,
                       CandDplusData const& candidatesDplus,
                       TracksWithExtra const& tracks)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus, tracks);
  }
  PROCESS_SWITCH(HfTaskDirectedFlowCharmHadrons, processDplusStd, "Process Dplus candidates with rectangular cuts", true);

  // Dstar with ML
  void processDstarMl(CollsWithQvecs::iterator const& collision,
                      CandDstarDataWMl const& candidatesDstar,
                      TracksWithExtra const& tracks)
  {
    runFlowAnalysis<DecayChannel::DstarToD0Pi>(collision, candidatesDstar, tracks);
  }
  PROCESS_SWITCH(HfTaskDirectedFlowCharmHadrons, processDstarMl, "Process Dstar candidates with ML", false);

  // Dstar with rectangular cuts
  void processDstarStd(CollsWithQvecs::iterator const& collision,
                       CandDstarData const& candidatesDstar,
                       TracksWithExtra const& tracks)
  {
    runFlowAnalysis<DecayChannel::DstarToD0Pi>(collision, candidatesDstar, tracks);
  }
  PROCESS_SWITCH(HfTaskDirectedFlowCharmHadrons, processDstarStd, "Process Dstar candidates with rectangular cuts", true);

}; // End struct HfTaskDirectedFlowCharmHadrons

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDirectedFlowCharmHadrons>(cfgc)};
}
