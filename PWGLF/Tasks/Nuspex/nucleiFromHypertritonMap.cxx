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
///
/// \author Roberta Ferioli (roberta.ferioli@cern.ch)
/// \since November, 2024

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TDatabasePDG.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TVector2.h>
#include <TVector3.h>

#include <string>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

namespace
{
constexpr int DeuteronPDG = o2::constants::physics::kDeuteron;
constexpr int HePDG = o2::constants::physics::kHelium3;
constexpr int HypertritonPDG = o2::constants::physics::kHyperTriton;
constexpr int HyperHelium4PDG = o2::constants::physics::kHyperHelium4;
} // namespace

struct nucleiFromHypertritonMap {

  using MCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using MCTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::McTrackLabels>;
  Preslice<aod::McParticles> mMcParticlesPerCollision = o2::aod::mcparticle::mcCollisionId;

  HistogramRegistry registryMC{
    "registryMC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Track Parameters
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 7, "minimum number of found ITS clusters"};
  Configurable<int> min_ITS_InnerBarrel_nClusters{"min_ITS_InnerBarrel_nClusters", 1, "minimum number of found ITS Inner Barrel clusters"};
  Configurable<int> max_ITS_InnerBarrel_nClusters{"max_ITS_InnerBarrel_nClusters", 3, "maximum number of found ITS Inner Barrel clusters"};
  Configurable<int> min_TPC_nClusters{"min_TPC_nClusters", 100, "minimum number of found TPC clusters"};
  Configurable<int> min_TPC_nCrossedRows{"min_TPC_nCrossedRows", 70, "minimum number of TPC crossed pad rows"};
  Configurable<float> max_chi2_TPC{"max_chi2_TPC", 4.0f, "maximum TPC chi^2/Ncls"};
  Configurable<float> min_chi2_TPC{"min_chi2_ITS", 0.5f, "minimum TPC chi^2/Ncls"};
  Configurable<float> min_eta{"min_eta", -0.8f, "minimum_eta"};
  Configurable<float> max_eta{"max_eta", +0.8f, "maximum_eta"};
  Configurable<float> max_dcaxy{"max_dcaxy", 0.05f, "Maximum DCAxy"};
  Configurable<float> max_dcaz{"max_dcaz", 0.05f, "Maximum DCAz"};
  Configurable<float> min_nsigmaTPC{"min_nsigmaTPC", -2.0f, "Minimum nsigma TPC"};
  Configurable<float> max_nsigmaTPC{"max_nsigmaTPC", +2.0f, "Maximum nsigma TPC"};

  Configurable<float> settingCutVertex{"settingCutVertex", 10.0f, "Accepted z-vertex range"};
  ConfigurableAxis pt_axis{"pt_axis", {50, -10.f, 10.f}, ";signed #it{p}_{T} (GeV/#it{c});Counts"};

  Configurable<bool> saveHelium{"saveHelium", false, "Save helium candidates"};

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgCCDBurl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> mCcdb;
  int mRunNumber = 0;
  float mBz = 0.f;

  int mSelectedPDG = 0;
  o2::dataformats::DCA mDcaInfoCov;
  o2::track::TrackParametrizationWithError<float> mTrackParCov;
  o2::base::MatLayerCylSet* mLut = nullptr;

  void init(InitContext const&)
  {
    mCcdb->setURL(cfgCCDBurl);
    mCcdb->setCaching(true);
    mCcdb->setLocalObjectValidityChecking();
    mCcdb->setFatalWhenNull(false);
    mLut = o2::base::MatLayerCylSet::rectifyPtrFromFile(mCcdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    registryMC.add("hypertritonPtGen", "hypertritonPtGen", HistType::kTH1F, {pt_axis});

    if (saveHelium) {
      registryMC.add("he3SecPtRec_from_hypertriton", "he3SecPtRec_from_hypertriton", HistType::kTH1F, {pt_axis});
      registryMC.add("he3SecPtGen_from_hypertriton", "he3SecPtGen_from_hypertriton", HistType::kTH1F, {pt_axis});
      registryMC.add("hyperHe4PtGen", "hyperHe4PtGen", HistType::kTH1F, {pt_axis});
      registryMC.add("he3SecPtRec_from_hyperHe4", "he3SecPtRec_from_hyperHe4", HistType::kTH1F, {pt_axis});
      registryMC.add("he3SecPtGen_from_hyperHe4", "he3SecPtGen_from_hyperHe4", HistType::kTH1F, {pt_axis});
      registryMC.add("he3PtRec", "he3PtRec", HistType::kTH1F, {pt_axis});
      registryMC.add("he3PtGen", "he3PtGen", HistType::kTH1F, {pt_axis});
    } else {
      registryMC.add("deutSecPtRec_from_hypertriton", "deutSecPtRec_from_hypertriton", HistType::kTH1F, {pt_axis});
      registryMC.add("deutSecPtGen_from_hypertriton", "deutSecPtGen_from_hypertriton", HistType::kTH1F, {pt_axis});
      registryMC.add("deutPtRec", "deutPtRec", HistType::kTH1F, {pt_axis});
      registryMC.add("deutPtGen", "deutPtGen", HistType::kTH1F, {pt_axis});
    }

    if (saveHelium) {
      mSelectedPDG = HePDG;
    } else {
      mSelectedPDG = DeuteronPDG;
    }
  }

  void initCCDB(const aod::BCsWithTimestamps::iterator& bc)
  {
    if (mRunNumber == bc.runNumber())
      return;

    auto timestamp = bc.timestamp();
    mRunNumber = bc.runNumber();

    o2::parameters::GRPMagField* grpmag = mCcdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(mLut);
    mBz = static_cast<float>(grpmag->getNominalL3Field());
    LOGF(info, "Retrieved GRP for timestamp %ull (%i) with magnetic field of %1.2f kZG", timestamp, mRunNumber, mBz);
  }

  void processMC(const MCCollisions& collisions, const aod::BCsWithTimestamps&, const aod::McParticles& mcParticles, const MCTracks& tracks)
  {
    for (const auto& collision : collisions) {

      if (!collision.sel8() || std::abs(collision.posZ()) > settingCutVertex) {
        return;
      }

      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      const int mcCollisionId = collision.mcCollisionId();
      auto mcParticlesThisCollision = mcParticles.sliceBy(mMcParticlesPerCollision, mcCollisionId);
      mcParticlesThisCollision.bindExternalIndices(&mcParticles);

      for (const auto& mcparticle : mcParticlesThisCollision) {

        if (((std::abs(mcparticle.pdgCode()) != HypertritonPDG && std::abs(mcparticle.pdgCode()) != HyperHelium4PDG) || !mcparticle.has_daughters()) && std::abs(mcparticle.pdgCode()) != mSelectedPDG) {
          continue;
        }

        const float particleSign = mcparticle.pdgCode() < 0 ? -1. : 1.;

        if (std::abs(mcparticle.pdgCode()) == HypertritonPDG) {

          for (const auto& daughter : mcparticle.daughters_as<aod::McParticles>()) {
            if (std::abs(daughter.pdgCode()) != mSelectedPDG) {
              continue;
            }
            const float daughterSign = daughter.pdgCode() < 0 ? -1. : 1.;

            registryMC.fill(HIST("hypertritonPtGen"), particleSign * mcparticle.pt());
            if (saveHelium) {
              registryMC.fill(HIST("he3SecPtGen_from_hypertriton"), daughterSign * daughter.pt());
            } else {
              registryMC.fill(HIST("deutSecPtGen_from_hypertriton"), daughterSign * daughter.pt());
            }
          }
        } else if (std::abs(mcparticle.pdgCode()) == HyperHelium4PDG) {
          for (const auto& daughter : mcparticle.daughters_as<aod::McParticles>()) {
            if (std::abs(daughter.pdgCode()) != mSelectedPDG) {
              continue;
            }
            const float daughterSign = daughter.pdgCode() < 0 ? -1. : 1.;

            registryMC.fill(HIST("hyperHe4PtGen"), particleSign * mcparticle.pt());
            if (saveHelium) {
              registryMC.fill(HIST("he3SecPtGen_from_hyperHe4"), daughterSign * daughter.pt());
            }
          }
        } else if (std::abs(mcparticle.pdgCode()) == HePDG && mcparticle.isPhysicalPrimary()) {
          registryMC.fill(HIST("he3PtGen"), particleSign * mcparticle.pt());
        } else if (std::abs(mcparticle.pdgCode()) == DeuteronPDG && mcparticle.isPhysicalPrimary()) {
          registryMC.fill(HIST("deutPtGen"), particleSign * mcparticle.pt());
        }
      }

      const o2::math_utils::Point3D<float> collisionVertex{collision.posX(), collision.posY(), collision.posZ()};

      for (const auto& track : tracks) {

        if (!track.has_mcParticle()) {
          continue;
        }

        const auto& mcparticle = track.mcParticle();
        if (std::abs(mcparticle.pdgCode()) != mSelectedPDG) {
          continue;
        }

        mDcaInfoCov.set(999, 999, 999, 999, 999);
        setTrackParCov(track, mTrackParCov);
        mTrackParCov.setPID(track.pidForTracking());
        std::array<float, 2> dcaInfo;
        o2::base::Propagator::Instance()->propagateToDCA(collisionVertex, mTrackParCov, mBz, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfo);

        const float dcaXY = dcaInfo[0];
        const float dcaZ = dcaInfo[1];

        if (track.itsNCls() < min_ITS_nClusters ||
            track.itsNClsInnerBarrel() < min_ITS_InnerBarrel_nClusters ||
            track.itsNClsInnerBarrel() > max_ITS_InnerBarrel_nClusters ||
            track.tpcNClsFound() < min_TPC_nClusters ||
            track.tpcNClsCrossedRows() < min_TPC_nCrossedRows ||
            track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
            track.tpcChi2NCl() > max_chi2_TPC ||
            track.tpcChi2NCl() < min_chi2_TPC ||
            track.eta() < min_eta || track.eta() > max_eta ||
            dcaXY > max_dcaxy || dcaXY < -max_dcaxy ||
            dcaZ > max_dcaz || dcaZ < -max_dcaz ||
            track.itsChi2NCl() > 36.f) {
          continue;
        }

        const float particleSign = mcparticle.pdgCode() < 0 ? -1. : 1.;

        if (std::abs(mcparticle.pdgCode()) == DeuteronPDG && mcparticle.isPhysicalPrimary()) {
          registryMC.fill(HIST("deutPtRec"), particleSign * track.pt());
        }
        if (std::abs(mcparticle.pdgCode()) == HePDG && mcparticle.isPhysicalPrimary()) {
          registryMC.fill(HIST("he3PtRec"), particleSign * 2 * track.pt());
        }

        for (const auto& motherparticle : mcparticle.mothers_as<aod::McParticles>()) {

          if (std::abs(motherparticle.pdgCode()) == HypertritonPDG) {
            if (std::abs(mcparticle.pdgCode()) == HePDG) {
              registryMC.fill(HIST("he3SecPtRec_from_hypertriton"), 2 * particleSign * track.pt());
            } else {
              registryMC.fill(HIST("deutSecPtRec_from_hypertriton"), particleSign * track.pt());
            }
          } else if (std::abs(motherparticle.pdgCode()) == HyperHelium4PDG) {
            if (std::abs(mcparticle.pdgCode()) == HePDG) {
              registryMC.fill(HIST("he3SecPtRec_from_hyperHe4"), 2 * particleSign * track.pt());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(nucleiFromHypertritonMap, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nucleiFromHypertritonMap>(cfgc)};
}
