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

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <vector>
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "PHOSBase/Geometry.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DetectorsBase/Propagator.h"

/// \struct PHOS electron id analysis
/// \brief Task for calculating electron identification parameters
/// \author Yeghishe Hambardzumyan, MIPT
/// \since Apr, 2024
/// @note Inherits functions and variables from phosAlign & phosPi0.
/// @note Results will be used for candidate table producing task.
///

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phosElId {

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using tracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;

  Configurable<float> mMinE{"mMinCluE", 0.3, "Minimum cluster energy for analysis"};

  Configurable<int> nBinsDeltaX{"nBinsDeltaX", 250, "N bins for track and cluster coordinate delta"};
  Configurable<float> mDeltaXmin{"mDeltaXmin", -50., "Min for track and cluster coordinate delta"};
  Configurable<float> mDeltaXmax{"mDeltaXmax", 50., "Max for track and cluster coordinate delta"};

  Configurable<int> nBinsDeltaZ{"nBinsDeltaZ", 250, "N bins for track and cluster coordinate delta"};
  Configurable<float> mDeltaZmin{"mDeltaZmin", -50., "Min for track and cluster coordinate delta"};
  Configurable<float> mDeltaZmax{"mDeltaZmax", 50., "Max for track and cluster coordinate delta"};

  Configurable<int> nBinsEp{"nBinsEp", 400, "N bins for E/p histograms"};
  Configurable<float> mEpmin{"mEpmin", -1., "Min for E/p histograms"};
  Configurable<float> mEpmax{"mEpmax", 3., "Max for E/p histograms"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::unique_ptr<o2::phos::Geometry> geomPHOS;
  double bz{0.}; // magnetic field
  int runNumber{0};

  HistogramRegistry mHistManager{"phosElIdHistograms"};

  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS electron identification analysis task ...";

    std::vector<double> momentum_binning = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0,
                                            1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                            4.5, 5.0, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.};
    const AxisSpec
      axisP{momentum_binning, "p (GeV/c)"},
      axisCounter{1, 0, +1, ""},
      axisEp{nBinsEp, mEpmin, mEpmax, "E/p", "E_{cluster}/p_{track}"},
      axisdX{nBinsDeltaX, mDeltaXmin, mDeltaXmax, "x_{tr}-x_{clu} (cm)", "x_{tr}-x_{clu} (cm)"},
      axisdZ{nBinsDeltaZ, mDeltaZmin, mDeltaZmax, "z_{tr}-z_{clu} (cm)", "z_{tr}-z_{clu} (cm)"},
      axisModes{4, 1., 5., "module", "module"},
      axisX{150, -75., 75., "x (cm)", "x (cm)"},
      axisZ{150, -75., 75., "z (cm)", "z (cm)"};

    mHistManager.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    mHistManager.add("hDeltaXPos_v_p", "positive trackX - clusterX vs p", HistType::kTH2F, {axisdX, axisP});
    mHistManager.add("hDeltaXNeg_v_p", "negative trackX - clusterX vs p", HistType::kTH2F, {axisdX, axisP});
    mHistManager.add("hDeltaZ_v_p", "trackZ - clusterZ vs p", HistType::kTH2F, {axisdZ, axisP});
    mHistManager.add("hEp_v_p", "E/p ratio vs p", HistType::kTH2F, {axisEp, axisP});
    mHistManager.add("hEp_v_p_disp", "E/p ratio vs p | OK dispersion", HistType::kTH2F, {axisEp, axisP});
    mHistManager.add("hdXdZp", "dx,dz,p_{tr}", HistType::kTH3F, {axisdX, axisdZ, axisP});
    mHistManager.add("hdZpmod", "dz,p_{tr},module", HistType::kTH3F, {axisdZ, axisP, axisModes});
    mHistManager.add("hdXpmod", "dx,p_{tr},module", HistType::kTH3F, {axisdX, axisP, axisModes});
    mHistManager.add("hdXpmod_pos", "dx,p_{tr},module positive tracks", HistType::kTH3F, {axisdX, axisP, axisModes});
    mHistManager.add("hdXpmod_neg", "dx,p_{tr},module negative tracks", HistType::kTH3F, {axisdX, axisP, axisModes});

    geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");
  }
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               aod::CaloClusters& clusters,
               tracks& tracks,
               aod::BCsWithTimestamps const&)
  {
    mHistManager.fill(HIST("eventCounter"), 0.5);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      runNumber = bc.runNumber();
    }

    if (clusters.size() == 0)
      return; // Nothing to process

    for (auto const& track : tracks) {
      if (!track.has_collision())
        continue;
      if (std::abs(track.collision_as<SelCollisions>().posZ()) > 10.f)
        continue;

      // calculate coordinate in PHOS plane
      if (std::abs(track.eta()) > 0.3)
        continue;
      int16_t module;
      float trackX = 999., trackZ = 999.;

      auto trackPar = getTrackPar(track);
      if (!impactOnPHOS(trackPar, track.trackEtaEmcal(), track.trackPhiEmcal(), track.collision_as<SelCollisions>().posZ(), module, trackX, trackZ))
        continue;

      float trackMom = track.p();
      for (const auto& clu : clusters) {
        if (clu.e() < mMinE)
          continue;
        if (module != clu.mod())
          continue;
        bool isDispOK = testLambda(clu.e(), clu.m02(), clu.m20());
        float posX = clu.x(), posZ = clu.z(), dX = trackX - posX, dZ = trackZ - posZ;

        if (track.sign() > 0)
          mHistManager.fill(HIST("hDeltaXPos_v_p"), dX, trackMom);
        else
          mHistManager.fill(HIST("hDeltaXNeg_v_p"), dX, trackMom);
        mHistManager.fill(HIST("hDeltaZ_v_p"), dZ, trackMom);
        mHistManager.fill(HIST("hEp_v_p"), clu.e() / trackMom, trackMom);
        mHistManager.fill(HIST("hdZpmod"), dZ, trackMom, module);
        mHistManager.fill(HIST("hdXpmod"), dX, trackMom, module);
        if (track.sign() > 0)
          mHistManager.fill(HIST("hdXpmod_pos"), dX, trackMom, module);
        else
          mHistManager.fill(HIST("hdXpmod_neg"), dX, trackMom, module);

        if (isDispOK)
          mHistManager.fill(HIST("hEp_v_p_disp"), clu.e() / trackMom, trackMom);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////
  bool impactOnPHOS(o2::track::TrackParametrization<float>& trackPar, float trackEta, float trackPhi, float /*zvtx*/, int16_t& module, float& trackX, float& trackZ)
  {
    // eta,phi was calculated at EMCAL radius.
    // Extrapolate to PHOS assuming zeroB and current vertex
    // Check if direction in PHOS acceptance+20cm and return phos module number and coordinates in PHOS module plane
    const float phiMin = 240. * 0.017453293; // degToRad
    const float phiMax = 323. * 0.017453293; // PHOS+20 cm * degToRad
    const float etaMax = 0.178266;
    if (trackPhi < phiMin || trackPhi > phiMax || abs(trackEta) > etaMax) {
      return false;
    }

    const float dphi = 20. * 0.017453293;
    if (trackPhi < 0.) {
      trackPhi += TMath::TwoPi();
    }
    if (trackPhi > TMath::TwoPi()) {
      trackPhi -= TMath::TwoPi();
    }
    module = 1 + static_cast<int16_t>((trackPhi - phiMin) / dphi);
    if (module < 1) {
      module = 1;
    }
    if (module > 4) {
      module = 4;
    }

    // get PHOS radius
    constexpr float shiftY = -1.26;    // Depth-optimized
    double posL[3] = {0., 0., shiftY}; // local position at the center of module
    double posG[3] = {0};
    geomPHOS->getAlignmentMatrix(module)->LocalToMaster(posL, posG);
    double rPHOS = sqrt(posG[0] * posG[0] + posG[1] * posG[1]);
    double alpha = (230. + 20. * module) * 0.017453293;

    // During main reconstruction track was propagated to radius 460 cm with accounting material
    // now material is not available. Therefore, start from main rec. position and extrapoate to actual radius without material
    // Get track parameters at point where main reconstruction stop
    float xPHOS = 460.f, xtrg = 0.f;
    if (!trackPar.getXatLabR(xPHOS, xtrg, bz, o2::track::DirType::DirOutward)) {
      return false;
    }
    auto prop = o2::base::Propagator::Instance();
    if (!trackPar.rotate(alpha) ||
        !prop->PropagateToXBxByBz(trackPar, xtrg, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrNONE)) {
      return false;
    }
    // calculate xyz from old (Phi, eta) and new r
    float r = std::sqrt(trackPar.getX() * trackPar.getX() + trackPar.getY() * trackPar.getY());
    trackPar.setX(r * std::cos(trackPhi - alpha));
    trackPar.setY(r * std::sin(trackPhi - alpha));
    trackPar.setZ(r / std::tan(2. * std::atan(std::exp(-trackEta))));

    if (!prop->PropagateToXBxByBz(trackPar, rPHOS, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrNONE)) {
      return false;
    }
    alpha = trackPar.getAlpha();
    double ca = cos(alpha), sa = sin(alpha);
    posG[0] = trackPar.getX() * ca - trackPar.getY() * sa;
    posG[1] = trackPar.getY() * ca + trackPar.getX() * sa;
    posG[2] = trackPar.getZ();

    geomPHOS->getAlignmentMatrix(module)->MasterToLocal(posG, posL);
    trackX = posL[0];
    trackZ = posL[1];
    return true;
  }
  //_____________________________________________________________________________
  bool testLambda(float pt, float l1, float l2)
  {
    // Parameterization for full dispersion
    // Parameterizatino for full dispersion
    float l2Mean = 1.53126 + 9.50835e+06 / (1. + 1.08728e+07 * pt + 1.73420e+06 * pt * pt);
    float l1Mean = 1.12365 + 0.123770 * TMath::Exp(-pt * 0.246551) + 5.30000e-03 * pt;
    float l2Sigma = 6.48260e-02 + 7.60261e+10 / (1. + 1.53012e+11 * pt + 5.01265e+05 * pt * pt) + 9.00000e-03 * pt;
    float l1Sigma = 4.44719e-04 + 6.99839e-01 / (1. + 1.22497e+00 * pt + 6.78604e-07 * pt * pt) + 9.00000e-03 * pt;
    float c = -0.35 - 0.550 * TMath::Exp(-0.390730 * pt);

    return 0.5 * (l1 - l1Mean) * (l1 - l1Mean) / l1Sigma / l1Sigma +
             0.5 * (l2 - l2Mean) * (l2 - l2Mean) / l2Sigma / l2Sigma +
             0.5 * c * (l1 - l1Mean) * (l2 - l2Mean) / l1Sigma / l2Sigma <
           4.;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{
    adaptAnalysisTask<phosElId>(cfgc)};
  return workflow;
}
