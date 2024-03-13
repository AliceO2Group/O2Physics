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
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/EventPlaneHelper.h"
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
  kFV0A,
  kBpos,
  kBneg,
  kNqVecDetectors
};
static const std::vector<std::string> qVecDetectorNames{"FT0C", "FT0A", "FV0A", "Bpos", "Bneg"};

std::shared_ptr<TH3> hQxQy[kNqVecDetectors];
std::shared_ptr<TH3> hNormQxQy[kNqVecDetectors];
std::shared_ptr<TH2> hPsi[kNqVecDetectors];
std::shared_ptr<TH2> hDeltaPsi[kNqVecDetectors][kNqVecDetectors];
std::shared_ptr<TH2> hScalarProduct[kNqVecDetectors][kNqVecDetectors];
std::shared_ptr<TH2> hNormalisedScalarProduct[kNqVecDetectors][kNqVecDetectors];
} // namespace

struct flowQC {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "Centrality estimator (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3)"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};

  ConfigurableAxis cfgCentralityBins{"cfgCentralityBins", {100, 0., 100.}, "Centrality binning"};
  ConfigurableAxis cfgQvecBins{"cfgQvecBins", {100, -2.f, 2.f}, "Binning for scalar product"};
  ConfigurableAxis cfgPhiBins{"cfgPhiBins", {140, -3.5f, 3.5f}, "Binning for azimuthal angle"};
  ConfigurableAxis cfgDeltaPhiBins{"cfgDeltaPhiBins", {280, -7.f, 7.f}, "Binning for azimuthal-angle differences"};

  // CCDB options
  Configurable<double> cfgBz{"cfgBz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> cfgCCDBurl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> cfgGRPpath{"cfgGRPpath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> cfgGRPmagPath{"cfgGRPmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  int mRunNumber = 0;
  float mBz = 0.f;

  // Flow analysis
  using CollWithQvec = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>::iterator;

  HistogramRegistry flow{"flow", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  template <class collision_t>
  bool eventSelection(collision_t& collision)
  {
    return collision.sel8() && collision.posZ() > -cfgCutVertex && collision.posZ() < cfgCutVertex;
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

    const AxisSpec psiAxis{cfgPhiBins, "#psi"};

    // z vertex histogram
    flow.add("hRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});

    // Centrality histograms
    flow.add("hCentFT0C", "", HistType::kTH1F, {centFT0cAxis});
    flow.add("hCentFT0A", "", HistType::kTH1F, {centFT0aAxis});
    flow.add("hCentFT0M", "", HistType::kTH1F, {centFT0mAxis});
    flow.add("hCentFV0A", "", HistType::kTH1F, {centFV0aAxis});

    for (int iQvecDet = 0; iQvecDet < qVecDetectors::kNqVecDetectors; iQvecDet++) {
      hQxQy[iQvecDet] = flow.add<TH3>(Form("hQxQy_%s", qVecDetectorNames[iQvecDet].c_str()), "", HistType::kTH3F, {centAxis, QxAxis, QyAxis});
      hNormQxQy[iQvecDet] = flow.add<TH3>(Form("hNormQxQy_%s", qVecDetectorNames[iQvecDet].c_str()), "", HistType::kTH3F, {centAxis, NormQxAxis, NormQyAxis});
      hPsi[iQvecDet] = flow.add<TH2>(Form("hPsi_%s", qVecDetectorNames[iQvecDet].c_str()), "", HistType::kTH2F, {centAxis, psiAxis});
      for (int jQvecDet = iQvecDet + 1; jQvecDet < qVecDetectors::kNqVecDetectors; jQvecDet++) {

        // Q-vector azimuthal-angle differences
        hDeltaPsi[iQvecDet][jQvecDet] = flow.add<TH2>(Form("hDeltaPsi_%s_%s", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str()), "", HistType::kTH2F, {centAxis, {cfgDeltaPhiBins, Form("#psi_{%s} - #psi_{%s}", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str())}});

        // Scalar-product histograms
        auto spLabel = Form("#vec{Q}_{2}^{%s} #upoint #vec{Q}_{2}^{%s}", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str());

        hScalarProduct[iQvecDet][jQvecDet] = flow.add<TH2>(Form("hScalarProduct_%s_%s", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str()), "", HistType::kTH2F, {centAxis, {cfgQvecBins, spLabel}});

        // Normalised scalar-product histograms
        auto normSpLabel = Form("#frac{#vec{Q}_{2}^{%s} #upoint #vec{Q}_{2}^{%s}}{||#vec{Q}_{2}^{%s}|| ||#vec{Q}_{2}^{%s}||}", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str(), qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str());

        hNormalisedScalarProduct[iQvecDet][jQvecDet] = flow.add<TH2>(Form("hNormalisedScalarProduct_%s_%s", qVecDetectorNames[iQvecDet].c_str(), qVecDetectorNames[jQvecDet].c_str()), "", HistType::kTH2F, {centAxis, {cfgQvecBins, normSpLabel}});
      }
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

  void process(CollWithQvec const& collision, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    gRandom->SetSeed(bc.timestamp());

    flow.fill(HIST("hRecVtxZData"), collision.posZ());

    const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};

    flow.fill(HIST("hCentFT0C"), collision.centFT0C());
    flow.fill(HIST("hCentFT0A"), collision.centFT0A());
    flow.fill(HIST("hCentFT0M"), collision.centFT0M());
    flow.fill(HIST("hCentFV0A"), collision.centFV0A());

    float centrality = getCentrality(collision);

    float QxFT0C = collision.qvecFT0CRe();
    float QyFT0C = collision.qvecFT0CIm();
    float QmodFT0C = std::hypot(QxFT0C, QyFT0C);
    float psiFT0C = std::atan2(QyFT0C, QxFT0C) / 2;

    float QxFT0A = collision.qvecFT0ARe();
    float QyFT0A = collision.qvecFT0AIm();
    float QmodFT0A = std::hypot(QxFT0A, QyFT0A);
    float psiFT0A = std::atan2(QyFT0A, QxFT0A) / 2;

    float QxFV0A = collision.qvecFV0ARe();
    float QyFV0A = collision.qvecFV0AIm();
    float QmodFV0A = std::hypot(QxFV0A, QyFV0A);
    float psiFV0A = std::atan2(QyFV0A, QxFV0A) / 2;

    float QxBpos = collision.qvecBPosRe();
    float QyBpos = collision.qvecBPosIm();
    float QmodBpos = std::hypot(QxBpos, QyBpos);
    float psiBpos = std::atan2(QyBpos, QxBpos) / 2;

    float QxBneg = collision.qvecBNegRe();
    float QyBneg = collision.qvecBNegIm();
    float QmodBneg = std::hypot(QxBneg, QyBneg);
    float psiBneg = std::atan2(QyBneg, QxBneg) / 2;

    std::array<float, qVecDetectors::kNqVecDetectors> vec_Qx = {QxFT0C, QxFT0A, QxFV0A, QxBpos, QxBneg};
    std::array<float, qVecDetectors::kNqVecDetectors> vec_Qy = {QyFT0C, QyFT0A, QyFV0A, QyBpos, QyBneg};
    std::array<float, qVecDetectors::kNqVecDetectors> vec_Qmod = {QmodFT0C, QmodFT0A, QmodFV0A, QmodBpos, QmodBneg};
    std::array<float, qVecDetectors::kNqVecDetectors> vec_Qpsi = {psiFT0C, psiFT0A, psiFV0A, psiBpos, psiBneg};

    for (int iQvecDet = 0; iQvecDet < qVecDetectors::kNqVecDetectors; iQvecDet++) {
      hQxQy[iQvecDet]->Fill(centrality, vec_Qx[iQvecDet], vec_Qy[iQvecDet]);
      hNormQxQy[iQvecDet]->Fill(centrality, vec_Qx[iQvecDet] / vec_Qmod[iQvecDet], vec_Qy[iQvecDet] / vec_Qmod[iQvecDet]);
      hPsi[iQvecDet]->Fill(centrality, vec_Qpsi[iQvecDet]);
      for (int jQvecDet = iQvecDet + 1; jQvecDet < qVecDetectors::kNqVecDetectors; jQvecDet++) {
        // Q-vector azimuthal-angle differences
        hDeltaPsi[iQvecDet][jQvecDet]->Fill(centrality, vec_Qpsi[iQvecDet] - vec_Qpsi[jQvecDet]);
        // Scalar-product histograms
        auto getSP = [&](int iDet1, int iDet2) {
          return vec_Qx[iDet1] * vec_Qx[iDet2] + vec_Qy[iDet1] * vec_Qy[iDet2];
        };
        hScalarProduct[iQvecDet][jQvecDet]->Fill(centrality, getSP(iQvecDet, jQvecDet));
        // Normalised scalar-product histograms
        auto getNormSP = [&](int iDet1, int iDet2) {
          return getSP(iDet1, iDet2) / (vec_Qmod[iDet1] * vec_Qmod[iDet2]);
        };
        hNormalisedScalarProduct[iQvecDet][jQvecDet]->Fill(centrality, getNormSP(iQvecDet, jQvecDet));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowQC>(cfgc, TaskName{"flow-qc"})};
}
