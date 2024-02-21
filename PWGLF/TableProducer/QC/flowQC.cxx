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
  kFV0A = 0,
  kFT0M = 1,
  kFT0A = 2,
  kFT0C = 3
};
static const std::vector<std::string> centDetectorNames{"FV0A", "FT0M", "FT0A", "FT0C"};
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

    const AxisSpec centAxis{cfgCentralityBins, fmt::format("{} percentile", (std::string)centDetectorNames[cfgCentralityEstimator])};

    const AxisSpec QxAxis{cfgQvecBins, "Q_{x}"};
    const AxisSpec QyAxis{cfgQvecBins, "Q_{y}"};

    const AxisSpec psiAxis{cfgPhiBins, "#psi"};

    const AxisSpec deltaPsiFT0cFT0a{cfgDeltaPhiBins, "#psi_{FT0C} - #psi_{FT0A}"};
    const AxisSpec deltaPsiFV0aFT0a{cfgDeltaPhiBins, "#psi_{FV0A} - #psi_{FT0A}"};
    const AxisSpec deltaPsiFV0aFT0c{cfgDeltaPhiBins, "#psi_{FV0A} - #psi_{FT0C}"};

    const AxisSpec ft0Aft0CspAxis{cfgQvecBins, "#vec{Q}_{2}^{FT0A} #upoint #vec{Q}_{2}^{FT0C}"};
    const AxisSpec fv0Aft0CspAxis{cfgQvecBins, "#vec{Q}_{2}^{FV0A} #upoint #vec{Q}_{2}^{FT0C}"};
    const AxisSpec fv0Aft0AspAxis{cfgQvecBins, "#vec{Q}_{2}^{FV0A} #upoint #vec{Q}_{2}^{FT0A}"};

    // z vertex histogram
    flow.add("hRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});

    // Centrality histograms
    flow.add("hCentFT0C", "", HistType::kTH1F, {centAxis});
    flow.add("hCentFT0A", "", HistType::kTH1F, {centAxis});
    flow.add("hCentFT0M", "", HistType::kTH1F, {centAxis});
    flow.add("hCentFV0A", "", HistType::kTH1F, {centAxis});

    // Q-vector histograms
    flow.add("hQxQyFT0C", "", HistType::kTH3F, {centAxis, QxAxis, QyAxis});
    flow.add("hQxQyFT0A", "", HistType::kTH3F, {centAxis, QxAxis, QyAxis});
    flow.add("hQxQyFV0A", "", HistType::kTH3F, {centAxis, QxAxis, QyAxis});

    // Q-vector azimuthal angles
    flow.add("hPsiFT0C", "", HistType::kTH2F, {centAxis, psiAxis});
    flow.add("hPsiFT0A", "", HistType::kTH2F, {centAxis, psiAxis});
    flow.add("hPsiFV0A", "", HistType::kTH2F, {centAxis, psiAxis});

    // Q-vector azimuthal-angle differences
    flow.add("hDeltaPsiFT0CFT0A", "", HistType::kTH2F, {centAxis, deltaPsiFT0cFT0a});
    flow.add("hDeltaPsiFV0AFT0A", "", HistType::kTH2F, {centAxis, deltaPsiFV0aFT0a});
    flow.add("hDeltaPsiFV0AFT0C", "", HistType::kTH2F, {centAxis, deltaPsiFV0aFT0c});

    // Scalar-product histograms
    flow.add("hScalarProductFT0AvsFT0C", "", HistType::kTH2F, {centAxis, ft0Aft0CspAxis});
    flow.add("hScalarProductFV0AvsFT0C", "", HistType::kTH2F, {centAxis, fv0Aft0CspAxis});
    flow.add("hScalarProductFV0AvsFT0A", "", HistType::kTH2F, {centAxis, fv0Aft0AspAxis});
  }

  template <typename Tcoll>
  float getCentrality(Tcoll const& collision)
  {
    if (cfgCentralityEstimator == centDetectors::kFV0A) {
      return collision.centFV0A();
    } else if (cfgCentralityEstimator == centDetectors::kFT0M) {
      return collision.centFT0M();
    } else if (cfgCentralityEstimator == centDetectors::kFT0A) {
      return collision.centFT0A();
    } else if (cfgCentralityEstimator == centDetectors::kFT0C) {
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
    float psiFT0C = std::atan2(QyFT0C, QxFT0C);

    float QxFT0A = collision.qvecFT0ARe();
    float QyFT0A = collision.qvecFT0AIm();
    float psiFT0A = std::atan2(QyFT0A, QxFT0A);

    float QxFV0A = collision.qvecFV0ARe();
    float QyFV0A = collision.qvecFV0AIm();
    float psiFV0A = std::atan2(QyFV0A, QxFV0A);

    flow.fill(HIST("hQxQyFT0C"), centrality, QxFT0C, QyFT0C);
    flow.fill(HIST("hQxQyFT0A"), centrality, QxFT0A, QyFT0A);
    flow.fill(HIST("hQxQyFV0A"), centrality, QxFV0A, QyFV0A);

    flow.fill(HIST("hPsiFT0C"), centrality, psiFT0C);
    flow.fill(HIST("hPsiFT0A"), centrality, psiFT0A);
    flow.fill(HIST("hPsiFV0A"), centrality, psiFV0A);

    flow.fill(HIST("hDeltaPsiFT0CFT0A"), centrality, psiFT0C - psiFT0A);
    flow.fill(HIST("hDeltaPsiFV0AFT0A"), centrality, psiFV0A - psiFT0A);
    flow.fill(HIST("hDeltaPsiFV0AFT0C"), centrality, psiFV0A - psiFT0C);

    flow.fill(HIST("hScalarProductFT0AvsFT0C"), centrality, QxFT0A * QxFT0C + QyFT0A * QyFT0C);
    flow.fill(HIST("hScalarProductFV0AvsFT0C"), centrality, QxFV0A * QxFT0C + QyFV0A * QyFT0C);
    flow.fill(HIST("hScalarProductFV0AvsFT0A"), centrality, QxFT0A * QxFV0A + QyFT0A * QyFV0A);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowQC>(cfgc, TaskName{"flow-qc"})};
}
