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

/// \file taskv1CharmHadrons.cxx
/// \brief Analysis task for charm hadron flow
///
/// \author Prottay Das, prottay.das@cern.ch

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "TF1.h"

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

enum DecayChannel { DplusToPiKPi = 0 };

struct HfTaskDirectedFlowCharmHadrons {
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<int> selectionFlag{"selectionFlag", 1, "Selection Flag for hadron (e.g. 1 for skimming, 3 for topo. and kine., 7 for PID)"};
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality accepted in SP computation"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality accepted in SP computation"};
  Configurable<bool> storeMl{"storeMl", false, "Flag to store ML scores"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indexes of BDT scores to be stored. Two indexes max."};

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
  using CollsWithQvecs = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::SPCalibrationTables>;
  using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;

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
    registry.add("hpuxQxpvscentpteta", "hpuxQxpvscentpteta", HistType::kTHnSparseF, axes, true);
    registry.add("hpuyQypvscentpteta", "hpuyQypvscentpteta", HistType::kTHnSparseF, axes, true);
    registry.add("hpuxQxtvscentpteta", "hpuxQxtvscentpteta", HistType::kTHnSparseF, axes, true);
    registry.add("hpuyQytvscentpteta", "hpuyQytvscentpteta", HistType::kTHnSparseF, axes, true);

    registry.add("hpQxtQxpvscent", "hpQxtQxpvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
    registry.add("hpQytQypvscent", "hpQytQypvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
    registry.add("hpQxtQypvscent", "hpQxtQypvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);
    registry.add("hpQxpQytvscent", "hpQxpQytvscent", HistType::kTHnSparseF, {thnAxisCent, thnAxisScalarProd}, true);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }; // end init

  /// Fill THnSparse
  /// \param mass is the invariant mass of the candidate
  /// \param pt is the transverse momentum of the candidate
  /// \param cent is the centrality of the collision
  /// \param cosNPhi is the cosine of the n*phi angle
  /// \param cosDeltaPhi is the cosine of the n*(phi - evtPl) angle
  /// \param sp is the scalar product
  /// \param outputMl are the ML scores
  void fillThn(double& mass,
               double& pt,
               double& eta,
               double& cent,
               double& uxQxp,
               double& uyQyp,
               double& uxQxt,
               double& uyQyt,
               std::vector<double>& outputMl,
               double& sign)
  {
    if (storeMl) {
      registry.fill(HIST("hpuxQxpvscentpteta"), mass, cent, pt, eta, uxQxp, sign, outputMl[0], outputMl[1]);
      registry.fill(HIST("hpuyQypvscentpteta"), mass, cent, pt, eta, uyQyp, sign, outputMl[0], outputMl[1]);
      registry.fill(HIST("hpuxQxtvscentpteta"), mass, cent, pt, eta, uxQxt, sign, outputMl[0], outputMl[1]);
      registry.fill(HIST("hpuyQytvscentpteta"), mass, cent, pt, eta, uyQyt, sign, outputMl[0], outputMl[1]);
    } else {
      registry.fill(HIST("hpuxQxpvscentpteta"), mass, cent, pt, eta, uxQxp, sign);
      registry.fill(HIST("hpuyQypvscentpteta"), mass, cent, pt, eta, uyQyp, sign);
      registry.fill(HIST("hpuxQxtvscentpteta"), mass, cent, pt, eta, uxQxt, sign);
      registry.fill(HIST("hpuyQytvscentpteta"), mass, cent, pt, eta, uyQyt, sign);
    }
  }

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

  /// Check if the collision is selected
  /// \param collision is the collision with the Q vector information
  /// \param bc is the bunch crossing with timestamp information
  /// \return true if the collision is selected, false otherwise
  template <o2::hf_centrality::CentralityEstimator centEstimator>
  bool isCollSelected(CollsWithQvecs::iterator const& collision,
                      aod::BCsWithTimestamps const&)
  {
    double centrality{-1.f};
    const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

    /// monitor the satisfied event selections
    hfEvSel.fillHistograms(collision, rejectionMask, centrality);
    return rejectionMask == 0;
  }

  /// Get the Q vector
  /// \param collision is the collision with the Q vector information
  std::vector<double> getQvec(CollsWithQvecs::iterator const& collision)
  {
    double qxAQVec = 0.0;
    double qyAQVec = 0.0;
    double qxCQVec = 0.0;
    double qyCQVec = 0.0;
    qxAQVec = collision.qxZDCA();
    qyAQVec = collision.qyZDCA();
    qxCQVec = collision.qxZDCC();
    qyCQVec = collision.qyZDCC();
    return {qxAQVec, qyAQVec, qxCQVec, qyCQVec};
  }

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param candidates are the selected candidates
  template <DecayChannel channel, typename T1, typename Trk>
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision,
                       T1 const& candidates, Trk const& /*tracks*/)
  {
    double cent = getCentrality(collision);
    if (cent < centralityMin || cent > centralityMax) {
      return;
    }

    if (!collision.triggerevent()) { // for selecting only callibrated events
      return;
    }

    std::vector<double> qVecs = getQvec(collision);
    double qxAQVec = qVecs[0];
    double qyAQVec = qVecs[1];
    double qxCQVec = qVecs[2];
    double qyCQVec = qVecs[3];

    auto qxZDCA = qxAQVec;
    auto qyZDCA = qyAQVec;
    auto qxZDCC = qxCQVec;
    auto qyZDCC = qyCQVec;

    auto QxtQxp = qxZDCC * qxZDCA;
    auto QytQyp = qyZDCC * qyZDCA;
    auto QxpQyt = qxZDCA * qyZDCC;
    auto QxtQyp = qxZDCC * qyZDCA;

    // denominators for SP calculation
    registry.fill(HIST("hpQxtQxpvscent"), cent, QxtQxp);
    registry.fill(HIST("hpQytQypvscent"), cent, QytQyp);
    registry.fill(HIST("hpQxpQytvscent"), cent, QxpQyt);
    registry.fill(HIST("hpQxtQypvscent"), cent, QxtQyp);

    for (const auto& candidate : candidates) {
      double massCand = 0.;
      std::vector<double> outputMl = {-999., -999.};
      if constexpr (std::is_same_v<T1, CandDplusData> || std::is_same_v<T1, CandDplusDataWMl>) {
        massCand = hfHelper.invMassDplusToPiKPi(candidate);
        if constexpr (std::is_same_v<T1, CandDplusDataWMl>) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++)
            outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
        }
      }

      auto trackprong0 = candidate.template prong0_as<Trk>();
      double sign = trackprong0.sign();

      double ptCand = candidate.pt();
      double etaCand = candidate.eta();
      double phiCand = candidate.phi();
      double cosNPhi = std::cos(phiCand);
      double sinNPhi = std::sin(phiCand);

      auto ux = cosNPhi;
      auto uy = sinNPhi;
      auto uxQxp = ux * qxZDCA;
      auto uyQyp = uy * qyZDCA;
      auto uxQxt = ux * qxZDCC;
      auto uyQyt = uy * qyZDCC;

      fillThn(massCand, ptCand, etaCand, cent, uxQxp, uyQyp, uxQxt, uyQyt, outputMl, sign);
    }
  }
  // Dplus with ML
  void processDplusMl(CollsWithQvecs::iterator const& collision,
                      CandDplusDataWMl const& candidatesDplus, TracksWithExtra const& tracks)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus, tracks);
  }
  PROCESS_SWITCH(HfTaskDirectedFlowCharmHadrons, processDplusMl, "Process Dplus candidates with ML", false);

  // Dplus with rectangular cuts
  void processDplus(CollsWithQvecs::iterator const& collision,
                    CandDplusData const& candidatesDplus, TracksWithExtra const& tracks)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus, tracks);
  }
  PROCESS_SWITCH(HfTaskDirectedFlowCharmHadrons, processDplus, "Process Dplus candidates", true);

}; // End struct HfTaskDirectedFlowCharmHadrons

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDirectedFlowCharmHadrons>(cfgc)};
}
