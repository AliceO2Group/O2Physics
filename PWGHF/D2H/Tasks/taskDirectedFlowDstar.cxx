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

/// \file taskDirectedFlowDstar.cxx
/// \brief Analysis task for Dstar directed flow
///
/// \author Sourav Kundu, sourav.kundudas@cern.ch

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

struct HfTaskDirectedFlowDstar {
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<bool> selectionFlag{"selectionFlag", true, "Selection Flag for D* decay to D0 Pi"};
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality accepted in SP computation"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality accepted in SP computation"};
  Configurable<bool> storeMl{"storeMl", false, "Flag to store ML scores"};
  Configurable<bool> direct{"direct", false, "Flag to calculate direct v1 odd and even"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indices of BDT scores to be stored. Two indexes max."};

  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {60, 0.138, 0.158}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0}, ""};
  ConfigurableAxis thnConfigAxisEta{"thnConfigAxisEta", {VARIABLE_WIDTH, -0.8, -0.4, 0, 0.4, 0.8}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0}, ""};
  ConfigurableAxis thnConfigAxisScalarProd{"thnConfigAxisScalarProd", {8000, -2.0, 2.0}, ""};
  ConfigurableAxis thnConfigAxisSign{"thnConfigAxisSign", {2, -2.0, 2.0}, ""};
  ConfigurableAxis thnConfigAxisMlOne{"thnConfigAxisMlOne", {1000, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisMlTwo{"thnConfigAxisMlTwo", {1000, 0., 1.}, ""};

  Filter filterSelectDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlag;

  using CandDstarDataWMl = soa::Filtered<soa::Join<aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi>>;
  using CandDstarData = soa::Filtered<soa::Join<aod::HfCandDstars, aod::HfSelDstarToD0Pi>>;
  using CollsWithQvecs = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::SPCalibrationTables>;
  using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

  SliceCache cache;
  HfHelper hfHelper;
  EventPlaneHelper epHelper;
  HfEventSelection hfEvSel; // event selection and monitoring
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
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
    } else {
      registry.add("hpQxtQxpvscent", "hpQxtQxpvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQytQypvscent", "hpQytQypvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQxtQypvscent", "hpQxtQypvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpQxpQytvscent", "hpQxpQytvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
      registry.add("hpuxQxpvscentpteta", "hpuxQxpvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpuyQypvscentpteta", "hpuyQypvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpuxQxtvscentpteta", "hpuxQxtvscentpteta", HistType::kTHnSparseF, axes, true);
      registry.add("hpuyQytvscentpteta", "hpuyQytvscentpteta", HistType::kTHnSparseF, axes, true);
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
  template <typename T1>
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision,
                       T1 const& candidates)
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
    } else {
      registry.fill(HIST("hpQxtQxpvscent"), cent, QxtQxp);
      registry.fill(HIST("hpQytQypvscent"), cent, QytQyp);
      registry.fill(HIST("hpQxpQytvscent"), cent, QxpQyt);
      registry.fill(HIST("hpQxtQypvscent"), cent, QxtQyp);
    }

    for (const auto& candidate : candidates) {
      double massCand = 0.;
      double sign = candidate.signSoftPi(); // soft pion sign
      std::vector<double> outputMl = {-999., -999.};
      if (sign > 0) {
        massCand = candidate.invMassDstar() - candidate.invMassD0();
      } else {
        massCand = candidate.invMassAntiDstar() - candidate.invMassD0Bar();
      }
      if constexpr (std::is_same_v<T1, CandDstarDataWMl>) {
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
          outputMl[iclass] = candidate.mlProbDstarToD0Pi()[classMl->at(iclass)];
      }
      double ptCand = candidate.pt();
      double etaCand = candidate.y(massCand);
      double phiCand = candidate.phi();
      double cosNPhi = std::cos(phiCand);
      double sinNPhi = std::sin(phiCand);

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
        } else {
          registry.fill(HIST("hpuxQxpvscentpteta"), massCand, cent, ptCand, etaCand, uxQxp, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpuyQypvscentpteta"), massCand, cent, ptCand, etaCand, uyQyp, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpuxQxtvscentpteta"), massCand, cent, ptCand, etaCand, uxQxt, sign, outputMl[0], outputMl[1]);
          registry.fill(HIST("hpuyQytvscentpteta"), massCand, cent, ptCand, etaCand, uyQyt, sign, outputMl[0], outputMl[1]);
        }
      } else {
        if (direct) {
          registry.fill(HIST("hpuxyQxypvscentpteta"), massCand, cent, ptCand, etaCand, uxyQxyp, sign);
          registry.fill(HIST("hpuxyQxytvscentpteta"), massCand, cent, ptCand, etaCand, uxyQxyt, sign);
          registry.fill(HIST("hpoddvscentpteta"), massCand, cent, ptCand, etaCand, oddv1, sign);
          registry.fill(HIST("hpevenvscentpteta"), massCand, cent, ptCand, etaCand, evenv1, sign);
        } else {
          registry.fill(HIST("hpuxQxpvscentpteta"), massCand, cent, ptCand, etaCand, uxQxp, sign);
          registry.fill(HIST("hpuyQypvscentpteta"), massCand, cent, ptCand, etaCand, uyQyp, sign);
          registry.fill(HIST("hpuxQxtvscentpteta"), massCand, cent, ptCand, etaCand, uxQxt, sign);
          registry.fill(HIST("hpuyQytvscentpteta"), massCand, cent, ptCand, etaCand, uyQyt, sign);
        }
      }
    }
  } // flow calculation
  // Dstar with ML
  void processDstarMl(CollsWithQvecs::iterator const& collision,
                      CandDstarDataWMl const& candidatesDstar)
  {
    runFlowAnalysis(collision, candidatesDstar);
  }
  PROCESS_SWITCH(HfTaskDirectedFlowDstar, processDstarMl, "Process Dstar candidates with ML", false);

  // Dstar with rectangular cuts
  void processDstarStd(CollsWithQvecs::iterator const& collision,
                       CandDstarData const& candidatesDstar)
  {
    runFlowAnalysis(collision, candidatesDstar);
  }
  PROCESS_SWITCH(HfTaskDirectedFlowDstar, processDstarStd, "Process Dstar candidates with rectangular cuts", true);

}; // End struct HfTaskDirectedFlowDstar

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDirectedFlowDstar>(cfgc)};
}
