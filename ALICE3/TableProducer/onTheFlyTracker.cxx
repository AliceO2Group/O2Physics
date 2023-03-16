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

/// \file onTheFlyTracker.cxx
///
/// \brief LUT-based on-the-fly analysis task-level tracking
///
/// This task allows for the calculation of aod::collisions and aod::Tracks in a synthetic manner,
/// smearing MC particles with very configurable settings. This will allow for the usage of
/// custom LUTs (obtained through separate studies) and the subsequent estimate of the performance
/// of a future detector even in very statistics-hungry analyses.
///
/// \author David Dobrigkeit Chinellato, UNICAMP

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"

using namespace o2;
using namespace o2::framework;

struct OnTheFlyTracker {
  Produces<aod::Collisions> collisions;
  Produces<aod::McCollisionLabels> collLabels;
  Produces<aod::StoredTracks> tracksPar;
  Produces<aod::TracksExtension> tracksParExtension;
  Produces<aod::StoredTracksCov> tracksParCov;
  Produces<aod::TracksCovExtension> tracksParCovExtension;
  Produces<aod::McTrackLabels> tracksLabels;
  Produces<aod::TracksDCA> tracksDCA;

  Configurable<float> maxEta{"maxEta", 1.5, "maximum eta to consider viable"};
  Configurable<float> minPt{"minPt", 0.1, "minimum pt to consider viable"};

  bool fillTracksDCA = false;

  // necessary for particle charges
  Service<O2DatabasePDG> pdg;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  void init(o2::framework::InitContext& initContext)
  {
    // Checking if the tables are requested in the workflow and enabling them
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      for (auto const& input : device.inputs) {
        if (input.matcher.binding == "TracksDCA") {
          fillTracksDCA = true;
        }
      }
    }

    // Basic QA
    const AxisSpec axisMomentum{static_cast<int>(100), 0.0f, +10.0f, "#it{p} (GeV/#it{c})"};
    histos.add("hPt", "hPt", kTH1F, {axisMomentum});
  }

  /// Function to convert a McParticle into a perfect Track
  /// \param particle the particle to convert (mcParticle)
  /// \param o2track the address of the resulting TrackParCov
  template <typename McParticleType>
  void convertMCParticleToO2Track(McParticleType& particle, o2::track::TrackParCov& o2track)
  {
    auto pdgInfo = pdg->GetParticle(particle.pdgCode());
    int charge = 0;
    if (pdgInfo != nullptr) {
      charge = pdgInfo->Charge();
    }
    std::array<float, 5> params;
    std::array<float, 15> covm = {0.};
    float s, c, x;
    o2::math_utils::sincos(particle.phi(), s, c);
    o2::math_utils::rotateZInv(particle.vx(), particle.vy(), x, params[0], s, c);
    params[1] = particle.vz();
    params[2] = 0.; // since alpha = phi
    auto theta = 2. * std::atan(std::exp(-particle.eta()));
    params[3] = 1. / std::tan(theta);
    params[4] = charge / particle.pt();

    // Initialize TrackParCov in-place
    new (&o2track)(o2::track::TrackParCov)(x, particle.phi(), params, covm);
  }

  /// Function to fill track parameter table
  /// \param coll collision (for index)
  /// \param trackType type of created track
  /// \param trackPar track for parameters
  template <typename CollType, typename TTrackPar>
  void fillTracksPar(CollType& coll, aod::track::TrackTypeEnum trackType, TTrackPar& trackPar)
  {
    tracksPar(coll.globalIndex(), trackType, trackPar.getX(), trackPar.getAlpha(), trackPar.getY(), trackPar.getZ(), trackPar.getSnp(), trackPar.getTgl(), trackPar.getQ2Pt());
    tracksParExtension(trackPar.getPt(), trackPar.getP(), trackPar.getEta(), trackPar.getPhi());
  }

  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    o2::dataformats::DCA dcaInfoCov;
    o2::dataformats::VertexBase vtx;

    for (const auto& mcParticle : mcParticles) {
      auto pdg = std::abs(mcParticle.pdgCode());
      if (pdg != kElectron && pdg != kMuonMinus && pdg != kPiPlus && pdg != kKPlus && pdg != kProton)
        continue;
      if (std::fabs(mcParticle.eta()) > maxEta)
        continue;
      if (mcParticle.pt() < minPt)
        continue;

      o2::track::TrackParCov trackParCov;
      convertMCParticleToO2Track(mcParticle, trackParCov);

      // *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
      // do the smearing in a one-liner with TrackSmearer
      // FIXME this has to be made available!
      // if (!smearer.smearTrack(o2track, pdg, nch)) continue;
      // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*

      // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
      // Calculate primary vertex
      // To be added once smeared tracks are in place
      // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*

      // Base QA
      histos.fill(HIST("hPt"), trackParCov.getPt());

      // Fixme: collision index could be changeable
      aod::track::TrackTypeEnum trackType = aod::track::Track;
      fillTracksPar(mcCollision, trackType, trackParCov);
      if (fillTracksDCA) {
        tracksDCA(1e-3, 1e-3);
      }
      // TODO do we keep the rho as 0? Also the sigma's are duplicated information
      tracksParCov(std::sqrt(trackParCov.getSigmaY2()), std::sqrt(trackParCov.getSigmaZ2()), std::sqrt(trackParCov.getSigmaSnp2()),
                   std::sqrt(trackParCov.getSigmaTgl2()), std::sqrt(trackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      tracksParCovExtension(trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(), trackParCov.getSigmaSnpY(),
                            trackParCov.getSigmaSnpZ(), trackParCov.getSigmaSnp2(), trackParCov.getSigmaTglY(), trackParCov.getSigmaTglZ(), trackParCov.getSigmaTglSnp(),
                            trackParCov.getSigmaTgl2(), trackParCov.getSigma1PtY(), trackParCov.getSigma1PtZ(), trackParCov.getSigma1PtSnp(), trackParCov.getSigma1PtTgl(),
                            trackParCov.getSigma1Pt2());
      tracksLabels(mcParticle.globalIndex(), 0);
    }
    collisions(-1, // BC is irrelevant in synthetic MC tests for now, could be adjusted in future
               mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
               1e-3, 0.0, 1e-3, 0.0, 0.0, 1e-3,
               0, 1e-3, mcParticles.size(),
               0, 0);
    collLabels(mcCollision.globalIndex(), 0);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyTracker>(cfgc)};
}
