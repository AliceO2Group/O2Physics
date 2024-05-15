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

#include <cmath>

#include "Math/Vector4D.h"

#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/EventPlaneHelper.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "Common/DataModel/Qvectors.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "TRandom3.h"

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
  kNqVecDetectors
};
static const std::vector<std::string> qVecDetectorNames{"FT0C", "FT0A", "TPCl", "TPCr"};

enum methods {
  kEP = 0,
  kQvec,
  kNmethods
};
static const std::vector<std::string> suffixes = {"EP", "Qvec"};

std::shared_ptr<TH3> hQxQy[kNmethods][kNqVecDetectors];
std::shared_ptr<TH3> hNormQxQy[kNmethods][kNqVecDetectors];
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

  // Flow analysis
  using CollWithEPandQvec = soa::Join<aod::Collisions,
                                      aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::EPCalibrationTables, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs>::iterator;

  HistogramRegistry general{"general", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry flow_ep{"flow_ep", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry flow_qvec{"flow_qvec", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry flow_comp{"flow_comp", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  template <class collision_t>
  bool eventSelection(collision_t& collision)
  {
    return collision.sel8() && collision.posZ() > -cfgCutVertex && collision.posZ() < cfgCutVertex && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.triggereventep() && collision.selection_bit(aod::evsel::kNoSameBunchPileup);
  }

  float computeEventPlane(float y, float x)
  {
    return 0.5 * TMath::ATan2(y, x);
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

    const AxisSpec QxAxis{cfgQvecBins, "Q_{2,x}"};
    const AxisSpec QyAxis{cfgQvecBins, "Q_{2,y}"};

    const AxisSpec NormQxAxis{cfgQvecBins, "#frac{Q_{2,x}}{||#vec{Q_{2}}||}"};
    const AxisSpec NormQyAxis{cfgQvecBins, "#frac{Q_{2,y}}{||#vec{Q_{2}}||}"};

    const AxisSpec psiAxis{cfgPhiBins, "#psi_{2}"};
    const AxisSpec psiCompAxis{cfgPhiBins, "#psi_{2}^{EP} - #psi_{2}^{Qvec}"};
    const AxisSpec cosPsiCompAxis{cfgCosPhiBins, "cos[2(#psi_{2}^{EP} - #psi_{2}^{Qvec})]"};

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
        hQxQy[iMethod][iQvecDet] = registry->add<TH3>(Form("hQxQy_%s_%s", qVecDetectorNames[iQvecDet].c_str(), suffixes[iMethod].c_str()), "", HistType::kTH3F, {centAxis, QxAxis, QyAxis});
        hNormQxQy[iMethod][iQvecDet] = registry->add<TH3>(Form("hNormQxQy_%s_%s", qVecDetectorNames[iQvecDet].c_str(), suffixes[iMethod].c_str()), "", HistType::kTH3F, {centAxis, NormQxAxis, NormQyAxis});
        hPsi[iMethod][iQvecDet] = registry->add<TH2>(Form("hPsi_%s_%s", qVecDetectorNames[iQvecDet].c_str(), suffixes[iMethod].c_str()), "", HistType::kTH2F, {centAxis, psiAxis});
        for (int jQvecDet = iQvecDet + 1; jQvecDet < qVecDetectors::kNqVecDetectors; jQvecDet++) {

          // Q-vector azimuthal-angle differences
          hDeltaPsi[iMethod][iQvecDet][jQvecDet] = registry->add<TH2>(Form("hDeltaPsi_%s_%s_%s", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str(), suffixes[iMethod].c_str()), "", HistType::kTH2F, {centAxis, {cfgDeltaPhiBins, Form("#psi_{%s} - #psi_{%s}", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str())}});

          // Scalar-product histograms
          auto spLabel = Form("#vec{Q}_{2}^{%s} #upoint #vec{Q}_{2}^{%s}", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str());

          hScalarProduct[iMethod][iQvecDet][jQvecDet] = registry->add<TH2>(Form("hScalarProduct_%s_%s_%s", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str(), suffixes[iMethod].c_str()), "", HistType::kTH2F, {centAxis, {cfgQvecBins, spLabel}});

          // Normalised scalar-product histograms
          auto normSpLabel = Form("#frac{#vec{Q}_{2}^{%s} #upoint #vec{Q}_{2}^{%s}}{||#vec{Q}_{2}^{%s}|| ||#vec{Q}_{2}^{%s}||}", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str(), qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str());

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

    // EP method
    float psiFT0A_EP = collision.psiFT0A();
    float QxFT0A_EP = std::cos(2 * psiFT0A_EP);
    float QyFT0A_EP = std::sin(2 * psiFT0A_EP);
    float QmodFT0A_EP = std::hypot(QxFT0A_EP, QyFT0A_EP);

    float psiFT0C_EP = collision.psiFT0C();
    float QxFT0C_EP = std::cos(2 * psiFT0C_EP);
    float QyFT0C_EP = std::sin(2 * psiFT0C_EP);
    float QmodFT0C_EP = std::hypot(QxFT0C_EP, QyFT0C_EP);

    float psiTPCl_EP = collision.psiTPCL();
    float QxTPCl_EP = std::cos(2 * psiTPCl_EP);
    float QyTPCl_EP = std::sin(2 * psiTPCl_EP);
    float QmodTPCl_EP = std::hypot(QxTPCl_EP, QyTPCl_EP);

    float psiTPCr_EP = collision.psiTPCR();
    float QxTPCr_EP = std::cos(2 * psiTPCr_EP);
    float QyTPCr_EP = std::sin(2 * psiTPCr_EP);
    float QmodTPCr_EP = std::hypot(QxTPCr_EP, QyTPCr_EP);

    // Qvec method
    float QxFT0A_Qvec = collision.qvecFT0ARe();
    float QyFT0A_Qvec = collision.qvecFT0AIm();
    float QmodFT0A_Qvec = std::hypot(QxFT0A_Qvec, QyFT0A_Qvec);
    float psiFT0A_Qvec = computeEventPlane(QyFT0A_Qvec, QxFT0A_Qvec);

    float QxFT0C_Qvec = collision.qvecFT0CRe();
    float QyFT0C_Qvec = collision.qvecFT0CIm();
    float QmodFT0C_Qvec = std::hypot(QxFT0C_Qvec, QyFT0C_Qvec);
    float psiFT0C_Qvec = computeEventPlane(QyFT0C_Qvec, QxFT0A_Qvec);

    float QxTPCl_Qvec = collision.qvecBNegRe();
    float QyTPCl_Qvec = collision.qvecBNegIm();
    float QmodTPCl_Qvec = std::hypot(QxTPCl_Qvec, QyTPCl_Qvec);
    float psiTPCl_Qvec = computeEventPlane(QyTPCl_Qvec, QxTPCl_Qvec);

    float QxTPCr_Qvec = collision.qvecBPosRe();
    float QyTPCr_Qvec = collision.qvecBPosIm();
    float QmodTPCr_Qvec = std::hypot(QxTPCr_Qvec, QyTPCr_Qvec);
    float psiTPCr_Qvec = computeEventPlane(QyTPCr_Qvec, QxTPCr_Qvec);

    std::array<float, qVecDetectors::kNqVecDetectors> vec_Qx[2] = {{QxFT0C_EP, QxFT0A_EP, QxTPCl_EP, QxTPCr_EP}, {QxFT0C_Qvec, QxFT0A_Qvec, QxTPCl_Qvec, QxTPCr_Qvec}};
    std::array<float, qVecDetectors::kNqVecDetectors> vec_Qy[2] = {{QyFT0C_EP, QyFT0A_EP, QyTPCl_EP, QyTPCr_EP}, {QyFT0C_Qvec, QyFT0A_Qvec, QyTPCl_Qvec, QyTPCr_Qvec}};
    std::array<float, qVecDetectors::kNqVecDetectors> vec_Qmod[2] = {{QmodFT0C_EP, QmodFT0A_EP, QmodTPCl_EP, QmodTPCr_EP}, {QmodFT0C_Qvec, QmodFT0A_Qvec, QmodTPCl_Qvec, QmodTPCr_Qvec}};
    std::array<float, qVecDetectors::kNqVecDetectors> vec_Qpsi[2] = {{psiFT0C_EP, psiFT0A_EP, psiTPCl_EP, psiTPCr_EP}, {psiFT0C_Qvec, psiFT0A_Qvec, psiTPCl_Qvec, psiTPCr_Qvec}};

    for (int iMethod = 0; iMethod < methods::kNmethods; iMethod++) {
      for (int iQvecDet = 0; iQvecDet < qVecDetectors::kNqVecDetectors; iQvecDet++) {
        hQxQy[iMethod][iQvecDet]->Fill(centrality, vec_Qx[iMethod][iQvecDet], vec_Qy[iMethod][iQvecDet]);
        hNormQxQy[iMethod][iQvecDet]->Fill(centrality, vec_Qx[iMethod][iQvecDet] / vec_Qmod[iMethod][iQvecDet], vec_Qy[iMethod][iQvecDet] / vec_Qmod[iMethod][iQvecDet]);
        hPsi[iMethod][iQvecDet]->Fill(centrality, vec_Qpsi[iMethod][iQvecDet]);
        for (int jQvecDet = iQvecDet + 1; jQvecDet < qVecDetectors::kNqVecDetectors; jQvecDet++) {
          // Q-vector azimuthal-angle differences
          hDeltaPsi[iMethod][iQvecDet][jQvecDet]->Fill(centrality, vec_Qpsi[iMethod][iQvecDet] - vec_Qpsi[iMethod][jQvecDet]);
          // Scalar-product histograms
          auto getSP = [&](int iDet1, int iDet2) {
            return vec_Qx[iMethod][iDet1] * vec_Qx[iMethod][iDet2] + vec_Qy[iMethod][iDet1] * vec_Qy[iMethod][iDet2];
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
      hCosPsiComp[iQvecDet]->Fill(centrality, std::cos(2 * (vec_Qpsi[methods::kEP][iQvecDet] - vec_Qpsi[methods::kQvec][iQvecDet])));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowQC>(cfgc, TaskName{"flow-qc"})};
}
