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
// Flow Qc analysis task
// ========================
//
// Executable + dependencies to run in local (+ maybe some conveters):
//
// Data (run3):
// o2-analysis-timestamp, o2-analysis-event-selection, o2-analysis-centrality-table,
// o2-analysis-multiplicity-table, o2-analysis-ft0-corrected-table, o2-analysis-track-propagation,
// o2-analysis-trackselection, o2-analysis-qvector-table, o2-analysis-lf-flow-qc

#include "PWGLF/DataModel/EPCalibrationTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TString.h>

#include <fmt/format.h>

#include <GPUROOTCartesianFwd.h>

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace
{
enum centDetectors {
  kCentFV0A = 0,
  kCentFT0M = 1,
  kCentFT0A = 2,
  kCentFT0C = 3
};
static const std::vector<std::string> centDetectorNames{"FV0A", "FT0M", "FT0A", "FT0C"};

enum qVecDetectors {
  kFT0C = 0,
  kFT0A,
  kTPCl,
  kTPCr,
  kTPC,
  kNqVecDetectors
};
static const std::vector<std::string> qVecDetectorNames{"FT0C", "FT0A", "TPCl", "TPCr", "TPC"};

enum methods {
  kEP = 0,
  kQvec,
  kNmethods
};
static const std::vector<std::string> suffixes = {"EP", "Qvec"};

std::shared_ptr<TH2> hPsi[kNmethods][kNqVecDetectors];
std::shared_ptr<TH2> hDeltaPsi[kNmethods][kNqVecDetectors][kNqVecDetectors];
std::shared_ptr<TH2> hScalarProduct[kNmethods][kNqVecDetectors][kNqVecDetectors];
std::shared_ptr<TH2> hNormalisedScalarProduct[kNmethods][kNqVecDetectors][kNqVecDetectors];

std::shared_ptr<TH2> hPsiComp[kNqVecDetectors];
std::shared_ptr<TH2> hCosPsiComp[kNqVecDetectors];
} // namespace

struct flowQC {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "Centrality estimator (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3)"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};

  ConfigurableAxis cfgCentralityBins{"cfgCentralityBins", {100, 0., 100.}, "Centrality binning"};
  ConfigurableAxis cfgQvecBins{"cfgQvecBins", {100, -2.f, 2.f}, "Binning for scalar product"};
  ConfigurableAxis cfgPhiBins{"cfgPhiBins", {140, -3.5f, 3.5f}, "Binning for azimuthal angle"};
  ConfigurableAxis cfgDeltaPhiBins{"cfgDeltaPhiBins", {280, -7.f, 7.f}, "Binning for azimuthal-angle differences"};
  ConfigurableAxis cfgCosPhiBins{"cfgCosPhiBins", {220, -1.1f, 1.1f}, "Binning for consinus of azimuthal angle"};

  // CCDB options
  Configurable<double> cfgBz{"cfgBz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> cfgCCDBurl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> cfgGRPpath{"cfgGRPpath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> cfgGRPmagPath{"cfgGRPmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  int mRunNumber = 0;
  float mBz = 0.f;

  Configurable<float> cfgHarmonic{"cfgHarmonic", 2.f, "Harmonics for flow analysis"};
  Configurable<bool> cfgQuadraticResponse{"cfgQuadraticResponse", false, "Use quadratic response for Q-vector quantities"};

  // Flow analysis
  using CollWithEPandQvec = soa::Join<aod::Collisions,
                                      aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::EPCalibrationTables, aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorFV0AVecs, aod::QvectorTPCallVecs, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs>::iterator;

  HistogramRegistry general{"general", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry flow_ep{"flow_ep", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry flow_qvec{"flow_qvec", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry flow_comp{"flow_comp", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  template <class collision_t>
  bool eventSelection(collision_t& collision)
  {
    return collision.sel8() && collision.posZ() > -cfgCutVertex && collision.posZ() < cfgCutVertex && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.triggereventep() && collision.selection_bit(aod::evsel::kNoSameBunchPileup);
  }

  float computeEventPlane(float y, float x, float harmonic)
  {
    return (1.f / harmonic) * TMath::ATan2(y, x);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();
    mRunNumber = bc.runNumber();

    if (cfgBz > -990) {
      mBz = cfgBz;
    } else {
      o2::parameters::GRPObject* grpo{ccdb->getForTimeStamp<o2::parameters::GRPObject>(cfgGRPpath, run3grp_timestamp)};
      o2::parameters::GRPMagField* grpmag{nullptr};
      if (grpo) {
        mBz = grpo->getNominalL3Field();
      } else {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(cfgGRPmagPath, run3grp_timestamp);
        if (!grpmag) {
          LOG(fatal) << "Got nullptr from CCDB for path " << cfgGRPmagPath << " of object GRPMagField and " << cfgGRPpath << " of object GRPObject for timestamp " << run3grp_timestamp;
        }
        mBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      }
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
    }
  }

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(cfgCCDBurl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    const AxisSpec centFT0aAxis{cfgCentralityBins, "FT0A percentile"};
    const AxisSpec centFT0cAxis{cfgCentralityBins, "FT0C percentile"};
    const AxisSpec centFT0mAxis{cfgCentralityBins, "FT0M percentile"};
    const AxisSpec centFV0aAxis{cfgCentralityBins, "FV0A percentile"};

    const AxisSpec centAxis{cfgCentralityBins, fmt::format("{} percentile", (std::string)centDetectorNames[cfgCentralityEstimator])};

    const char* qLabel = cfgQuadraticResponse ? "Q^{2}" : "Q";

    const AxisSpec psiAxis{cfgPhiBins, Form("#psi_{%.0f}", cfgHarmonic.value)};
    const AxisSpec psiCompAxis{cfgPhiBins, Form("#psi_{%.0f}^{EP} - #psi_{%.0f}^{Qvec}", cfgHarmonic.value, cfgHarmonic.value)};
    const AxisSpec cosPsiCompAxis{cfgCosPhiBins, Form("cos[2(#psi_{%.0f}^{EP} - #psi_{%.0f}^{Qvec})]", cfgHarmonic.value, cfgHarmonic.value)};

    // z vertex histogram
    general.add("hRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});

    // Centrality histograms
    general.add("hCentFT0C", "", HistType::kTH1F, {centFT0cAxis});
    general.add("hCentFT0A", "", HistType::kTH1F, {centFT0aAxis});
    general.add("hCentFT0M", "", HistType::kTH1F, {centFT0mAxis});
    general.add("hCentFV0A", "", HistType::kTH1F, {centFV0aAxis});

    for (int iMethod = 0; iMethod < methods::kNmethods; iMethod++) {
      HistogramRegistry* registry = (iMethod == methods::kEP) ? &flow_ep : &flow_qvec;

      for (int iQvecDet = 0; iQvecDet < qVecDetectors::kNqVecDetectors; iQvecDet++) {
        hPsi[iMethod][iQvecDet] = registry->add<TH2>(Form("hPsi_%s_%s", qVecDetectorNames[iQvecDet].c_str(), suffixes[iMethod].c_str()), "", HistType::kTH2F, {centAxis, psiAxis});
        for (int jQvecDet = iQvecDet + 1; jQvecDet < qVecDetectors::kNqVecDetectors; jQvecDet++) {

          // Q-vector azimuthal-angle differences
          hDeltaPsi[iMethod][iQvecDet][jQvecDet] = registry->add<TH2>(Form("hDeltaPsi_%s_%s_%s", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str(), suffixes[iMethod].c_str()), "", HistType::kTH2F, {centAxis, {cfgDeltaPhiBins, Form("#psi_{%s} - #psi_{%s}", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str())}});

          // Scalar-product histograms
          auto spLabel = Form("#vec{%s}_{%.0f}^{%s} #upoint #vec{%s}_{%.0f}^{%s}", qLabel, cfgHarmonic.value, qVecDetectorNames[iQvecDet].c_str(), qLabel, cfgHarmonic.value, qVecDetectorNames[jQvecDet].c_str());

          hScalarProduct[iMethod][iQvecDet][jQvecDet] = registry->add<TH2>(Form("hScalarProduct_%s_%s_%s", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str(), suffixes[iMethod].c_str()), "", HistType::kTH2F, {centAxis, {cfgQvecBins, spLabel}});

          // Normalised scalar-product histograms
          auto normSpLabel = Form("#frac{#vec{%s}_{%.0f}^{%s} #upoint #vec{%s}_{%.0f}^{%s}}{||#vec{%s}_{%.0f}^{%s}|| ||#vec{%s}_{%.0f}^{%s}||}", qLabel, cfgHarmonic.value, qVecDetectorNames[iQvecDet].c_str(), qLabel, cfgHarmonic.value, qVecDetectorNames[jQvecDet].c_str(), qLabel, cfgHarmonic.value, qVecDetectorNames[iQvecDet].c_str(), qLabel, cfgHarmonic.value, qVecDetectorNames[jQvecDet].c_str());

          hNormalisedScalarProduct[iMethod][iQvecDet][jQvecDet] = registry->add<TH2>(Form("hNormalisedScalarProduct_%s_%s_%s", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str(), suffixes[iMethod].c_str()), "", HistType::kTH2F, {centAxis, {cfgQvecBins, normSpLabel}});
        }
      }
    }
    for (int iQvecDet = 0; iQvecDet < qVecDetectors::kNqVecDetectors; iQvecDet++) {
      hPsiComp[iQvecDet] = flow_comp.add<TH2>(Form("hPsiComp_%s", qVecDetectorNames[iQvecDet].c_str()), "", HistType::kTH2F, {centAxis, psiCompAxis});
      hCosPsiComp[iQvecDet] = flow_comp.add<TH2>(Form("hCosPsiComp_%s", qVecDetectorNames[iQvecDet].c_str()), "", HistType::kTH2F, {centAxis, cosPsiCompAxis});
    }
  }

  template <typename Tcoll>
  float getCentrality(Tcoll const& collision)
  {
    if (cfgCentralityEstimator == centDetectors::kCentFV0A) {
      return collision.centFV0A();
    } else if (cfgCentralityEstimator == centDetectors::kCentFT0M) {
      return collision.centFT0M();
    } else if (cfgCentralityEstimator == centDetectors::kCentFT0A) {
      return collision.centFT0A();
    } else if (cfgCentralityEstimator == centDetectors::kCentFT0C) {
      return collision.centFT0C();
    } else {
      LOG(warning) << "Centrality estimator not valid. Possible values: (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3). Centrality set to 1.";
      return 1.;
    }
  }

  void process(CollWithEPandQvec const& collision, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    gRandom->SetSeed(bc.timestamp());

    if (!eventSelection(collision)) {
      return;
    }

    general.fill(HIST("hRecVtxZData"), collision.posZ());

    const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};

    general.fill(HIST("hCentFT0C"), collision.centFT0C());
    general.fill(HIST("hCentFT0A"), collision.centFT0A());
    general.fill(HIST("hCentFT0M"), collision.centFT0M());
    general.fill(HIST("hCentFV0A"), collision.centFV0A());

    float centrality = getCentrality(collision);

    auto maybeSquare = [this](float value) {
      return cfgQuadraticResponse ? value * value : value;
    };
    const float qvecHarmonic = cfgQuadraticResponse ? (cfgHarmonic.value / 2.f) : cfgHarmonic.value;
    const int qvecHarmonicIndex = static_cast<int>(qvecHarmonic) - 2;

    // EP method
    float QmodFT0A_EP = maybeSquare(collision.qFT0A());
    float psiFT0A_EP = collision.psiFT0A();

    float QmodFT0C_EP = maybeSquare(collision.qFT0C());
    float psiFT0C_EP = collision.psiFT0C();

    float QmodTPCl_EP = maybeSquare(collision.qTPCL());
    float psiTPCl_EP = collision.psiTPCL();

    float QmodTPCr_EP = maybeSquare(collision.qTPCR());
    float psiTPCr_EP = collision.psiTPCR();

    float QmodTPC_EP = maybeSquare(collision.qTPC());
    float psiTPC_EP = collision.psiTPC();

    // Qvec method
    float QxFT0A_Qvec_raw = collision.qvecFT0AReVec()[qvecHarmonicIndex];
    float QyFT0A_Qvec_raw = collision.qvecFT0AImVec()[qvecHarmonicIndex];
    float psiFT0A_Qvec = computeEventPlane(QyFT0A_Qvec_raw, QxFT0A_Qvec_raw, qvecHarmonic);
    float QmodFT0A_Qvec = maybeSquare(std::hypot(QxFT0A_Qvec_raw, QyFT0A_Qvec_raw));

    float QxFT0C_Qvec_raw = collision.qvecFT0CReVec()[qvecHarmonicIndex];
    float QyFT0C_Qvec_raw = collision.qvecFT0CImVec()[qvecHarmonicIndex];
    float psiFT0C_Qvec = computeEventPlane(QyFT0C_Qvec_raw, QxFT0C_Qvec_raw, qvecHarmonic);
    float QmodFT0C_Qvec = maybeSquare(std::hypot(QxFT0C_Qvec_raw, QyFT0C_Qvec_raw));

    float QxTPCl_Qvec_raw = collision.qvecTPCnegReVec()[qvecHarmonicIndex];
    float QyTPCl_Qvec_raw = collision.qvecTPCnegImVec()[qvecHarmonicIndex];
    float psiTPCl_Qvec = computeEventPlane(QyTPCl_Qvec_raw, QxTPCl_Qvec_raw, qvecHarmonic);
    float QmodTPCl_Qvec = maybeSquare(std::hypot(QxTPCl_Qvec_raw, QyTPCl_Qvec_raw));

    float QxTPCr_Qvec_raw = collision.qvecTPCposReVec()[qvecHarmonicIndex];
    float QyTPCr_Qvec_raw = collision.qvecTPCposImVec()[qvecHarmonicIndex];
    float psiTPCr_Qvec = computeEventPlane(QyTPCr_Qvec_raw, QxTPCr_Qvec_raw, qvecHarmonic);
    float QmodTPCr_Qvec = maybeSquare(std::hypot(QxTPCr_Qvec_raw, QyTPCr_Qvec_raw));

    float QxTPC_Qvec_raw = collision.qvecTPCallReVec()[qvecHarmonicIndex];
    float QyTPC_Qvec_raw = collision.qvecTPCallImVec()[qvecHarmonicIndex];
    float psiTPC_Qvec = computeEventPlane(QyTPC_Qvec_raw, QxTPC_Qvec_raw, qvecHarmonic);
    float QmodTPC_Qvec = maybeSquare(std::hypot(QxTPC_Qvec_raw, QyTPC_Qvec_raw));

    std::array<float, qVecDetectors::kNqVecDetectors> vec_Qmod[2] = {{QmodFT0C_EP, QmodFT0A_EP, QmodTPCl_EP, QmodTPCr_EP, QmodTPC_EP}, {QmodFT0C_Qvec, QmodFT0A_Qvec, QmodTPCl_Qvec, QmodTPCr_Qvec, QmodTPC_Qvec}};
    std::array<float, qVecDetectors::kNqVecDetectors> vec_Qpsi[2] = {{psiFT0C_EP, psiFT0A_EP, psiTPCl_EP, psiTPCr_EP, psiTPC_EP}, {psiFT0C_Qvec, psiFT0A_Qvec, psiTPCl_Qvec, psiTPCr_Qvec, psiTPC_Qvec}};

    for (int iMethod = 0; iMethod < methods::kNmethods; iMethod++) {
      for (int iQvecDet = 0; iQvecDet < qVecDetectors::kNqVecDetectors; iQvecDet++) {
        hPsi[iMethod][iQvecDet]->Fill(centrality, vec_Qpsi[iMethod][iQvecDet]);
        for (int jQvecDet = iQvecDet + 1; jQvecDet < qVecDetectors::kNqVecDetectors; jQvecDet++) {
          // Q-vector azimuthal-angle differences
          hDeltaPsi[iMethod][iQvecDet][jQvecDet]->Fill(centrality, vec_Qpsi[iMethod][iQvecDet] - vec_Qpsi[iMethod][jQvecDet]);
          // Scalar-product histograms
          auto getSP = [&](int iDet1, int iDet2) {
            return vec_Qmod[iMethod][iDet1] * vec_Qmod[iMethod][iDet2] * std::cos(cfgHarmonic.value * (vec_Qpsi[iMethod][iDet1] - vec_Qpsi[iMethod][iDet2]));
          };
          hScalarProduct[iMethod][iQvecDet][jQvecDet]->Fill(centrality, getSP(iQvecDet, jQvecDet));
          // Normalised scalar-product histograms
          auto getNormSP = [&](int iDet1, int iDet2) {
            return getSP(iDet1, iDet2) / (vec_Qmod[iMethod][iDet1] * vec_Qmod[iMethod][iDet2]);
          };
          hNormalisedScalarProduct[iMethod][iQvecDet][jQvecDet]->Fill(centrality, getNormSP(iQvecDet, jQvecDet));
        }
      }
    }
    for (int iQvecDet = 0; iQvecDet < qVecDetectors::kNqVecDetectors; iQvecDet++) {
      hPsiComp[iQvecDet]->Fill(centrality, vec_Qpsi[methods::kEP][iQvecDet] - vec_Qpsi[methods::kQvec][iQvecDet]);
      hCosPsiComp[iQvecDet]->Fill(centrality, std::cos(cfgHarmonic.value * (vec_Qpsi[methods::kEP][iQvecDet] - vec_Qpsi[methods::kQvec][iQvecDet])));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowQC>(cfgc, TaskName{"flow-qc"})};
}
