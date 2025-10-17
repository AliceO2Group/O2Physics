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
// Nuclei spectra analysis task
// ========================
//
// Executable + dependencies:
//
// Data (run3):
// o2-analysis-lf-nuclei-spectra, o2-analysis-timestamp
// o2-analysis-pid-tof-base, o2-analysis-multiplicity-table, o2-analysis-event-selection
// (to add flow: o2-analysis-qvector-table, o2-analysis-centrality-table)

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFSlimNucleiTables.h"
#include "PWGLF/Utils/nucleiUtils.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/Tools/TrackTuner.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector4D.h"
#include "TRandom3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct nucleiQC {

  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TOFSignal, aod::TOFEvTime, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCDe, aod::pidTOFDe, aod::pidTPCTr, aod::pidTOFTr, aod::pidTPCHe, aod::pidTOFHe, aod::pidTPCAl, aod::pidTOFAl>;
  using TrackCandidatesMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TOFSignal, aod::TOFEvTime, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCDe, aod::pidTOFDe, aod::pidTPCTr, aod::pidTOFTr, aod::pidTPCHe, aod::pidTOFHe, aod::pidTPCAl, aod::pidTOFAl, aod::McTrackLabels>;
  using Collision = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs, aod::McCollisionLabels>::iterator;
  Preslice<aod::McParticles> mMcParticlesPerCollision = o2::aod::mcparticle::mcCollisionId;

  Configurable<bool> cfgFillTable{"cfgFillTable", true, "Fill output tree"};
  Configurable<bool> cfgDoCheckPdgCode{"cfgDoCheckPdgCode", true, "Should you only select tracks associated to a mc particle with the correct PDG code?"};
  Configurable<LabeledArray<int>> cfgSpeciesToProcess{"cfgSpeciesToProcess", {nuclei::speciesToProcessDefault[0], nuclei::Species::kNspecies, 1, nuclei::names, {"processNucleus"}}, "Nuclei to process"};
  Configurable<LabeledArray<int>> cfgEventSelections{"cfgEventSelections", {nuclei::EvSelDefault[0], 8, 1, nuclei::eventSelectionLabels, nuclei::eventSelectionTitle}, "Event selections"};
  Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "Centrality estimator (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3)"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {nuclei::betheBlochDefault[0], nuclei::Species::kNspecies, 6, nuclei::names, nuclei::betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<int>> cfgUseCentralTpcCalibration{"cfgUseCentralTpcCalibration", {nuclei::useCentralTpcCalibrationDefault[0], nuclei::Species::kNspecies, 1, nuclei::names, {"UseCentralTpcCalibration"}}, "Use central TPC calibration"};

  Configurable<float> cfgRapidityMin{"cfgRapidityMin", -1., "Minimum rapidity value"};
  Configurable<float> cfgRapidityMax{"cfgRapidityMax", 1., "Maximum rapidity value"};
  Configurable<float> cfgRapidityCenterMass{"cfgRapidityCenterMass", 0.0f, "Center of mass rapidity"};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutTpcMom{"cfgCutTpcMom", 0.2f, "Minimum TPC momentum for tracks"};
  Configurable<float> cfgCutNclusITS{"cfgCutNclusITS", 5, "Minimum number of ITS clusters"};
  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};

  Configurable<LabeledArray<double>> cfgNsigmaTPC{"cfgNsigmaTPC", {nuclei::nSigmaTPCdefault[0], nuclei::Species::kNspecies, 2, nuclei::names, nuclei::nSigmaConfigName}, "TPC nsigma selection for light nuclei"};
  Configurable<LabeledArray<double>> cfgNsigmaTOF{"cfgNsigmaTOF", {nuclei::nSigmaTOFdefault[0], nuclei::Species::kNspecies, 2, nuclei::names, nuclei::nSigmaConfigName}, "TPC nsigma selection for light nuclei"};

  // CCDB options
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgCCDBurl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> mCcdb;
  int mRunNumber = 0;
  float mBz = 0.f;

  HistogramRegistry mHistograms{
    "histos",
    {
      {"hEventSelections", "Event selections; Selection step; Counts", {HistType::kTH1D, {{nuclei::evSel::kNevSels + 1, -0.5f, static_cast<float>(nuclei::evSel::kNevSels) + 0.5f}}}},
      {"hVtxZBefore", "Vertex distribution in Z before selections;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
      {"hVtxZ", "Vertex distribution in Z;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};
  std::vector<int> mSpeciesToProcess;
  Produces<aod::NucleiTableRed> mNucleiTableRed;

  std::vector<nuclei::SlimCandidate> mNucleiCandidates;
  std::vector<int> mFilledMcParticleIds;

  o2::dataformats::DCA mDcaInfoCov;
  o2::track::TrackParametrizationWithError<float> mTrackParCov;
  std::array<nuclei::PidManager, static_cast<int>(nuclei::Species::kNspecies)> mPidManagers;

  void init(o2::framework::InitContext&)
  {

    mCcdb->setURL(cfgCCDBurl);
    mCcdb->setCaching(true);
    mCcdb->setLocalObjectValidityChecking();
    mCcdb->setFatalWhenNull(false);
    nuclei::lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(mCcdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    for (int iSel = 0; iSel < nuclei::evSel::kNevSels; iSel++) {
      mHistograms.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(iSel + 1, nuclei::eventSelectionLabels[iSel].c_str());
    }

    for (int iSpecies = 0; iSpecies < static_cast<int>(nuclei::Species::kNspecies); iSpecies++) {
      if (cfgSpeciesToProcess->get(iSpecies) == 1) {
        mSpeciesToProcess.emplace_back(iSpecies);
      }
    }

    static_for<0, nuclei::kNspecies - 1>([&](auto iSpecies) {
      constexpr int kSpeciesCt = decltype(iSpecies)::value;
      const int kSpeciesRt = kSpeciesCt;

      if (std::find(mSpeciesToProcess.begin(), mSpeciesToProcess.end(), kSpeciesCt) == mSpeciesToProcess.end()) {
        return;
      }

      float tpcBetheBlochParams[6];
      for (int iParam = 0; iParam < 6; iParam++) {
        tpcBetheBlochParams[iParam] = cfgBetheBlochParams->get(kSpeciesRt, iParam);
      }

      nuclei::createHistogramRegistryNucleus<kSpeciesCt>(mHistograms);

      if (cfgUseCentralTpcCalibration->get(static_cast<uint32_t>(kSpeciesRt), static_cast<uint32_t>(0)) == 0) {
        mPidManagers[kSpeciesRt] = nuclei::PidManager(kSpeciesRt, tpcBetheBlochParams);
      } else {
        mPidManagers[kSpeciesRt] = nuclei::PidManager(kSpeciesRt);
      }
    });
  }

  void initCCDB(const aod::BCsWithTimestamps::iterator& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto timestamp = bc.timestamp();
    mRunNumber = bc.runNumber();

    o2::parameters::GRPMagField* grpmag = mCcdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(nuclei::lut);
    mBz = static_cast<float>(grpmag->getNominalL3Field());
    LOGF(info, "Retrieved GRP for timestamp %ull (%i) with magnetic field of %1.2f kZG", timestamp, mRunNumber, mBz);
  }

  template <typename Ttrack>
  bool trackSelection(const Ttrack& track)
  {
    if (std::abs(track.eta()) > cfgCutEta ||
        track.tpcInnerParam() < cfgCutTpcMom ||
        track.itsNCls() < cfgCutNclusITS ||
        track.tpcNClsFound() < cfgCutNclusTPC ||
        track.tpcNClsCrossedRows() < 70 ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcChi2NCl() > 4.f ||
        track.itsChi2NCl() > 36.f) {
      return false;
    }
    return true;
  }

  template <int iSpecies, typename Ttrack, typename Tcollision>
  bool pidSelection(const Ttrack& track, const Tcollision& collision)
  {
    constexpr int kIndex = iSpecies;
    if (!nuclei::checkSpeciesValidity(kIndex)) {
      std::runtime_error("species contains invalid nucleus kIndex");
    }

    float centrality = nuclei::getCentrality(collision, cfgCentralityEstimator);
    float nsigmaTPC = mPidManagers[kIndex].getNSigmaTPC(track);
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaTPC_preselectionVsCentrality"), track.pt() * track.sign(), nsigmaTPC, centrality);
    if (std::abs(nsigmaTPC) > cfgNsigmaTPC->get(kIndex, 1))
      return false;
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaTPCVsCentrality"), track.pt() * track.sign(), nsigmaTPC, centrality);

    float nsigmaITS = mPidManagers[kIndex].getNSigmaITS(track);
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaITS_preselectionVsCentrality"), track.sign() * track.pt(), nsigmaITS, centrality);
    // add nsigmaITS cut ?
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaITSVsCentrality"), track.sign() * track.pt(), nsigmaITS, centrality);

    float nsigmaTOF = mPidManagers[kIndex].getNSigmaTOF(track);
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaTOF_preselectionVsCentrality"), track.sign() * track.pt(), nsigmaTOF, centrality);
    if (std::abs(nsigmaTOF) > cfgNsigmaTOF->get(kIndex, 1) && track.hasTOF())
      return false;
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaTOFVsCentrality"), track.sign() * track.pt(), nsigmaTOF, centrality);

    return true;
  }

  template <typename Tparticle>
  void fillNucleusFlagsPdgsMc(const Tparticle& particle, nuclei::SlimCandidate& candidate)
  {
    candidate.pdgCode = particle.pdgCode();

    if (particle.isPhysicalPrimary()) {
      candidate.flags |= nuclei::Flags::kIsPhysicalPrimary;

      // heavy flavour mother
      // if (particle.has_mothers()) {
      //  for (const auto& motherparticle : particle.mothers_as<aod::McParticles>()) {
      //    if (std::find(nuclei::hfMothCodes.begin(), nuclei::hfMothCodes.end(), std::abs(motherparticle.pdgCode())) != nuclei::hfMothCodes.end()) {
      //      flags |= kIsSecondaryFromWeakDecay;
      //      motherPdgCode = motherparticle.pdgCode();
      //      break;
      //    }
      //  }
      //}

    } else if (particle.has_mothers()) {
      candidate.flags |= nuclei::Flags::kIsSecondaryFromWeakDecay;
      for (const auto& motherparticle : particle.template mothers_as<aod::McParticles>()) {
        candidate.motherPdgCode = motherparticle.pdgCode();
      }

    } else {
      candidate.flags |= nuclei::Flags::kIsSecondaryFromMaterial;
    }

    mFilledMcParticleIds.emplace_back(particle.globalIndex());
  }

  template <typename Tcollision, typename Ttrack>
  void fillNucleusFlagsPdgs(const int iSpecies, const Tcollision& collision, const Ttrack& track, nuclei::SlimCandidate& candidate)
  {
    candidate.flags = static_cast<uint16_t>((track.pidForTracking() & 0xF) << 12);
    candidate.flags |= iSpecies == nuclei::Species::kPr ? nuclei::Flags::kProton : iSpecies == nuclei::Species::kDe ? nuclei::Flags::kDeuteron
                                                                                 : iSpecies == nuclei::Species::kTr ? nuclei::Flags::kTriton
                                                                                 : iSpecies == nuclei::Species::kHe ? nuclei::Flags::kHe3
                                                                                 : iSpecies == nuclei::Species::kAl ? nuclei::Flags::kHe4
                                                                                                                    : 0;

    if (track.hasTOF()) {
      candidate.flags |= nuclei::Flags::kHasTOF;
    }
    if (track.hasTRD()) {
      candidate.flags |= nuclei::Flags::kHasTRD;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      candidate.flags |= nuclei::Flags::kITSrof;
    }
  }

  template <typename Tparticle>
  void fillNucleusGeneratedVariables(const Tparticle& particle, nuclei::SlimCandidate& candidate)
  {
    candidate.ptGenerated = particle.pt() * (particle.pdgCode() > 0 ? 1.f : -1.f);
    candidate.etaGenerated = particle.eta();
    candidate.phiGenerated = particle.phi();
  }

  template <typename Tcollision, typename Ttrack>
  void fillDcaInformation(const Tcollision& collision, const Ttrack& track, nuclei::SlimCandidate& candidate)
  {

    const o2::math_utils::Point3D<float> collisionVertex{collision.posX(), collision.posY(), collision.posZ()};

    mDcaInfoCov.set(999, 999, 999, 999, 999);
    setTrackParCov(track, mTrackParCov);
    mTrackParCov.setPID(track.pidForTracking());
    std::array<float, 2> dcaInfo;
    o2::base::Propagator::Instance()->propagateToDCA(collisionVertex, mTrackParCov, mBz, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfo);

    candidate.DCAxy = dcaInfo[0];
    candidate.DCAz = dcaInfo[1];
  }

  template <const bool isMc, typename Tcollision, typename Ttrack>
  nuclei::SlimCandidate fillCandidate(const int iSpecies, Tcollision const& collision, Ttrack const& track)
  {
    if (!nuclei::checkSpeciesValidity(iSpecies)) {
      std::runtime_error("species contains invalid nucleus index");
    }
    nuclei::SlimCandidate candidate = {.pt = track.pt() * track.sign(),
                                       .eta = track.eta(),
                                       .phi = track.phi(),
                                       .tpcInnerParam = track.tpcInnerParam(),
                                       .clusterSizesITS = track.itsClusterSizes(),
                                       .TPCsignal = track.tpcSignal(),
                                       .beta = mPidManagers[iSpecies].getBetaTOF(track),
                                       .DCAxy = 0.f,
                                       .DCAz = 0.f,
                                       .flags = 0,
                                       .pdgCode = 0,
                                       .motherPdgCode = 0,
                                       .ptGenerated = 0.f, // to be filled for mc
                                       .etaGenerated = 0.f,
                                       .phiGenerated = 0.f,
                                       .centrality = nuclei::getCentrality(collision, cfgCentralityEstimator)};

    fillDcaInformation(collision, track, candidate);
    fillNucleusFlagsPdgs(iSpecies, collision, track, candidate);

    if constexpr (isMc) {
      if (track.has_mcParticle()) {

        const auto& particle = track.mcParticle();
        fillNucleusFlagsPdgsMc(particle, candidate);
        fillNucleusGeneratedVariables(particle, candidate);
      }
    }

    mNucleiCandidates.emplace_back(candidate);
    return candidate;
  }

  template <bool isGenerated>
  void dispatchFillHistograms(const int iSpecies, const nuclei::SlimCandidate& candidate)
  {
    switch (iSpecies) {
      case nuclei::Species::kPr:
        return fillHistograms<nuclei::Species::kPr, isGenerated>(candidate);
      case nuclei::Species::kDe:
        return fillHistograms<nuclei::Species::kDe, isGenerated>(candidate);
      case nuclei::Species::kTr:
        return fillHistograms<nuclei::Species::kTr, isGenerated>(candidate);
      case nuclei::Species::kHe:
        return fillHistograms<nuclei::Species::kHe, isGenerated>(candidate);
      case nuclei::Species::kAl:
        return fillHistograms<nuclei::Species::kAl, isGenerated>(candidate);
      default:
        return;
    }
  }

  template <int iSpecies, const bool isGenerated>
  void fillHistograms(const nuclei::SlimCandidate& candidate)
  {
    constexpr int kIndex = iSpecies;
    if (!nuclei::checkSpeciesValidity(kIndex)) {
      std::runtime_error("species contains invalid nucleus kIndex");
    }

    if (isGenerated) {
      const float ptGenerated = (kIndex == nuclei::Species::kPr || kIndex == nuclei::Species::kDe || kIndex == nuclei::Species::kTr) ? candidate.ptGenerated : candidate.ptGenerated / 2.f;
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/hPtGenerated"), ptGenerated);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3PtVsEtaVsCentralityGenerated"), ptGenerated, candidate.etaGenerated, candidate.centrality);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3PhiVsEtaVsCentralityGenerated"), candidate.phiGenerated, candidate.etaGenerated, candidate.centrality);
    } else {
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/hPtReconstructed"), candidate.pt);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3PtVsEtaVsCentralityReconstructed"), candidate.pt, candidate.eta, candidate.centrality);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3PhiVsEtaVsCentralityReconstructed"), candidate.phi, candidate.eta, candidate.centrality);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3DCAxyVsPtVsCentrality"), candidate.pt, candidate.DCAxy, candidate.centrality);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3DCAzVsPtVsCentrality"), candidate.pt, candidate.DCAz, candidate.centrality);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3BetaVsPtVsCentrality"), candidate.pt, candidate.beta, candidate.centrality);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3dEdxVsPVsCentrality"), candidate.pt, candidate.TPCsignal, candidate.centrality);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3ClusterSizeVsPtVsCentrality"), candidate.pt, mPidManagers[kIndex].getClusterSizeCosLambdaITS(candidate.clusterSizesITS, candidate.eta), candidate.centrality);
    }
  }

  void processMc(const Collision& collision, const TrackCandidatesMC& tracks, const aod::BCsWithTimestamps&, const aod::McParticles& mcParticles)
  {
    mNucleiCandidates.clear();
    mFilledMcParticleIds.clear();

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!nuclei::eventSelection(collision, mHistograms, cfgEventSelections, cfgCutVertex)) {
      return;
    }

    for (const auto& track : tracks) {

      static_for<0, nuclei::kNspecies - 1>([&](auto iSpecies) {
        constexpr int kSpeciesCt = decltype(iSpecies)::value;
        const int kSpeciesRt = kSpeciesCt;

        if (std::find(mSpeciesToProcess.begin(), mSpeciesToProcess.end(), kSpeciesRt) == mSpeciesToProcess.end()) {
          return;
        }

        if (track.has_mcParticle()) {
          const auto& particle = track.mcParticle();
          if (cfgDoCheckPdgCode) {
            if (std::abs(particle.pdgCode()) != nuclei::pdgCodes[kSpeciesRt])
              return;
          }
          if ((particle.y() - cfgRapidityCenterMass) < cfgRapidityMin || (particle.y() - cfgRapidityCenterMass) > cfgRapidityMax) {
            return;
          }
        }

        mHistograms.fill(HIST(nuclei::cNames[kSpeciesCt]) + HIST("/hTrackSelections"), nuclei::trackSelection::kNoCuts);
        if (!trackSelection(track))
          return;
        mHistograms.fill(HIST(nuclei::cNames[kSpeciesCt]) + HIST("/hTrackSelections"), nuclei::trackSelection::kTrackCuts);

        if (!pidSelection<kSpeciesRt>(track, collision))
          return;
        mHistograms.fill(HIST(nuclei::cNames[kSpeciesCt]) + HIST("/hTrackSelections"), nuclei::trackSelection::kPidCuts);

        nuclei::SlimCandidate candidate;
        if (track.has_mcParticle()) {
          candidate = fillCandidate</*isMc*/ true>(kSpeciesCt, collision, track);
          dispatchFillHistograms</*isGenerated*/ true>(kSpeciesRt, candidate);
        } else {
          candidate = fillCandidate</*isMc*/ true>(kSpeciesCt, collision, track);
        }

        dispatchFillHistograms</*isGenerated*/ false>(kSpeciesRt, candidate);
      });
    }

    const int mcCollisionId = collision.mcCollisionId();
    auto mcParticlesThisCollision = mcParticles.sliceBy(mMcParticlesPerCollision, mcCollisionId);
    mcParticlesThisCollision.bindExternalIndices(&mcParticles);
    for (const auto& particle : mcParticlesThisCollision) {
      if (std::find(mFilledMcParticleIds.begin(), mFilledMcParticleIds.end(), particle.globalIndex()) != mFilledMcParticleIds.end()) {
        continue;
      }
      int iSpecies = nuclei::getSpeciesFromPdg(particle.pdgCode());
      if (std::find(mSpeciesToProcess.begin(), mSpeciesToProcess.end(), iSpecies) == mSpeciesToProcess.end()) {
        continue;
      }

      nuclei::SlimCandidate candidate;
      fillNucleusFlagsPdgsMc(particle, candidate);
      fillNucleusGeneratedVariables(particle, candidate);
      mNucleiCandidates.emplace_back(candidate);

      dispatchFillHistograms</*isGenerated*/ true>(iSpecies, candidate);
    }

    if (!cfgFillTable)
      return;

    for (const auto& candidate : mNucleiCandidates) {
      mNucleiTableRed(
        candidate.pt,
        candidate.eta,
        candidate.phi,
        candidate.tpcInnerParam,
        candidate.clusterSizesITS,
        candidate.TPCsignal,
        candidate.beta,
        candidate.DCAxy,
        candidate.DCAz,
        candidate.flags,
        candidate.pdgCode,
        candidate.motherPdgCode);
    }
  }
  PROCESS_SWITCH(nucleiQC, processMc, "Mc analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<nucleiQC>(cfgc)};
}
