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
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>, UniBo
/// \author Roberto Preghenella preghenella@bo.infn.it
///

#include <utility>

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
#include "TableHelper.h"

#include "ALICE3/Core/DelphesO2TrackSmearer.h"

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
  Configurable<float> multEtaRange{"multEtaRange", 0.8, "eta range to compute the multiplicity"};
  Configurable<float> minPt{"minPt", 0.1, "minimum pt to consider viable"};
  Configurable<bool> enableLUT{"enableLUT", false, "Enable track smearing"};
  Configurable<bool> enableNucleiSmearing{"enableNucleiSmearing", false, "Enable smearing of nuclei"};

  Configurable<std::string> lutEl{"lutEl", "lutCovm.el.dat", "LUT for electrons"};
  Configurable<std::string> lutMu{"lutMu", "lutCovm.mu.dat", "LUT for muons"};
  Configurable<std::string> lutPi{"lutPi", "lutCovm.pi.dat", "LUT for pions"};
  Configurable<std::string> lutKa{"lutKa", "lutCovm.ka.dat", "LUT for kaons"};
  Configurable<std::string> lutPr{"lutPr", "lutCovm.pr.dat", "LUT for protons"};
  Configurable<std::string> lutDe{"lutDe", "lutCovm.de.dat", "LUT for deuterons"};
  Configurable<std::string> lutTr{"lutTr", "lutCovm.tr.dat", "LUT for tritons"};
  Configurable<std::string> lutHe3{"lutHe3", "lutCovm.he3.dat", "LUT for Helium-3"};

  bool fillTracksDCA = false;

  // necessary for particle charges
  Service<O2DatabasePDG> pdgDB;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Track smearer
  o2::delphes::DelphesO2TrackSmearer mSmearer;

  void init(o2::framework::InitContext& initContext)
  {
    // Checking if the tables are requested in the workflow and enabling them
    fillTracksDCA = isTableRequiredInWorkflow(initContext, "TracksDCA");

    if (enableLUT) {
      std::map<int, const char*> mapPdgLut;
      //       const char* lutElChar = ((std::string)lutEl).c_str();
      //       const char* lutMuChar = ((std::string)lutMu).c_str();
      //       const char* lutPiChar = ((std::string)lutPi).c_str();
      //       const char* lutKaChar = ((std::string)lutKa).c_str();
      //       const char* lutPrChar = ((std::string)lutPr).c_str();
      mapPdgLut.insert(std::make_pair(11, "lutCovm.el.dat"));
      mapPdgLut.insert(std::make_pair(13, "lutCovm.mu.dat"));
      mapPdgLut.insert(std::make_pair(211, "lutCovm.pi.dat"));
      mapPdgLut.insert(std::make_pair(321, "lutCovm.ka.dat"));
      mapPdgLut.insert(std::make_pair(2212, "lutCovm.pr.dat"));
      if (enableNucleiSmearing) {
        const char* lutDeChar = ((std::string)lutDe).c_str();
        const char* lutTrChar = ((std::string)lutTr).c_str();
        const char* lutHe3Char = ((std::string)lutHe3).c_str();
        mapPdgLut.insert(std::make_pair(1000010020, lutDeChar));
        mapPdgLut.insert(std::make_pair(1000010030, lutTrChar));
        mapPdgLut.insert(std::make_pair(1000020030, lutHe3Char));
      }
      for (auto e : mapPdgLut) {
        if (!mSmearer.loadTable(e.first, e.second)) {
          LOG(fatal) << "Having issue with loading the LUT " << e.first << " " << e.second;
        }
      }
    }

    // Basic QA
    const AxisSpec axisMomentum{static_cast<int>(1000), 0.0f, +10.0f, "#it{p} (GeV/#it{c})"};
    histos.add("hPtGenerated", "hPtGenerated", kTH1F, {axisMomentum});
    histos.add("hPtGeneratedPi", "hPtGeneratedPi", kTH1F, {axisMomentum});
    histos.add("hPtGeneratedKa", "hPtGeneratedKa", kTH1F, {axisMomentum});
    histos.add("hPtGeneratedPr", "hPtGeneratedPr", kTH1F, {axisMomentum});
    histos.add("hPtReconstructed", "hPtReconstructed", kTH1F, {axisMomentum});
    histos.add("hPtReconstructedPi", "hPtReconstructedPi", kTH1F, {axisMomentum});
    histos.add("hPtReconstructedKa", "hPtReconstructedKa", kTH1F, {axisMomentum});
    histos.add("hPtReconstructedPr", "hPtReconstructedPr", kTH1F, {axisMomentum});
  }

  /// Function to convert a McParticle into a perfect Track
  /// \param particle the particle to convert (mcParticle)
  /// \param o2track the address of the resulting TrackParCov
  template <typename McParticleType>
  void convertMCParticleToO2Track(McParticleType& particle, o2::track::TrackParCov& o2track)
  {
    auto pdgInfo = pdgDB->GetParticle(particle.pdgCode());
    int charge = 0;
    if (pdgInfo != nullptr) {
      charge = pdgInfo->Charge() / 3;
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

  float dNdEta = 0.f; // Charged particle multiplicity to use in the efficiency evaluation
  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    o2::dataformats::DCA dcaInfoCov;
    o2::dataformats::VertexBase vtx;
    // First we compute the number of charged particles in the event
    dNdEta = 0.f;
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.eta()) > multEtaRange) {
        continue;
      }
      if (mcParticle.has_daughters()) {
        continue;
      }
      const auto& pdgInfo = pdgDB->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        LOG(warning) << "PDG code " << mcParticle.pdgCode() << " not found in the database";
        continue;
      }
      if (pdgInfo->Charge() == 0) {
        continue;
      }
      dNdEta += 1.f;
    }

    for (const auto& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      const auto pdg = std::abs(mcParticle.pdgCode());
      if (pdg != kElectron && pdg != kMuonMinus && pdg != kPiPlus && pdg != kKPlus && pdg != kProton) {
        continue;
      }
      if (std::fabs(mcParticle.eta()) > maxEta) {
        continue;
      }

      histos.fill(HIST("hPtGenerated"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 211)
        histos.fill(HIST("hPtGeneratedPi"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 321)
        histos.fill(HIST("hPtGeneratedKa"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 2212)
        histos.fill(HIST("hPtGeneratedPr"), mcParticle.pt());

      if (mcParticle.pt() < minPt) {
        continue;
      }
      o2::track::TrackParCov trackParCov;
      convertMCParticleToO2Track(mcParticle, trackParCov);

      if (!mSmearer.smearTrack(trackParCov, mcParticle.pdgCode(), dNdEta)) {
        continue;
      }

      // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
      // Calculate primary vertex
      // To be added once smeared tracks are in place
      // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*

      // Base QA (note: reco pT here)
      histos.fill(HIST("hPtReconstructed"), trackParCov.getPt());
      if (TMath::Abs(mcParticle.pdgCode()) == 211)
        histos.fill(HIST("hPtReconstructedPi"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 321)
        histos.fill(HIST("hPtReconstructedKa"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 2212)
        histos.fill(HIST("hPtReconstructedPr"), mcParticle.pt());

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

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<OnTheFlyTracker>(cfgc)}; }
