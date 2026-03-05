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
/// \brief TableProducer/Task for nuclei QC. The produced table can be disabled with a configurable.
///
/// \author Giorgio Alberto Lucia (giorgio.alberto.lucia@cern.ch)
///

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFSlimNucleiTables.h"
#include "PWGLF/Utils/nucleiUtils.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/Tools/TrackTuner.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/BetheBlochAleph.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector4D.h"
#include "TMCProcess.h"
#include "TRandom3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace
{

enum trackQuality {
  kNoCuts = 0,
  kEtaCut = 1,
  kPtpcCut = 2,
  kNclsItsCut = 3,
  kNclsTpcCut = 4,
  kNCrossedRowsCut = 5,
  kTpcChi2Cut = 6,
  kItsChi2Cut = 7,
  kNtrackQuality = 8
};

std::array<std::string, trackQuality::kNtrackQuality> trackQualityLabels{"All", "#eta cut", "#it{p}_{TPC}^{min} cut", "#it{N}_{cls}^{ITS} cut", "#it{N}_{cls}^{TPC} cut", "Crossed rows cut", "#chi^{2}_{TPC} cut", "#chi^{2}_{ITS} cut"};

} // namespace

struct nucleiQC {

  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TOFSignal, aod::TOFEvTime, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCDe, aod::pidTOFDe, aod::pidTPCTr, aod::pidTOFTr, aod::pidTPCHe, aod::pidTOFHe, aod::pidTPCAl, aod::pidTOFAl>;
  using TrackCandidatesMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TOFSignal, aod::TOFEvTime, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCDe, aod::pidTOFDe, aod::pidTPCTr, aod::pidTOFTr, aod::pidTPCHe, aod::pidTOFHe, aod::pidTPCAl, aod::pidTOFAl, aod::McTrackLabels>;
  using Collision = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs, aod::McCollisionLabels>::iterator;
  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs, aod::McCollisionLabels>;
  Preslice<TrackCandidatesMC> mTracksPerCollision = aod::track::collisionId;
  Preslice<aod::McParticles> mMcParticlesPerCollision = o2::aod::mcparticle::mcCollisionId;

  Configurable<bool> cfgFillTable{"cfgFillTable", true, "Fill output tree"};
  Configurable<bool> cfgDoCheckPdgCode{"cfgDoCheckPdgCode", true, "Should you only select tracks associated to a mc particle with the correct PDG code?"};
  Configurable<bool> cfgFillOnlyPhysicalPrimaries{"cfgFillOnlyPhysicalPrimaries", true, "Should you only select physical primary particles?"};
  Configurable<LabeledArray<int>> cfgSpeciesToProcess{"cfgSpeciesToProcess", {nuclei::speciesToProcessDefault[0], nuclei::Species::kNspecies, 1, nuclei::names, {"processNucleus"}}, "Nuclei to process"};
  Configurable<LabeledArray<int>> cfgEventSelections{"cfgEventSelections", {nuclei::EvSelDefault[0], 8, 1, nuclei::eventSelectionLabels, nuclei::eventSelectionTitle}, "Event selections"};
  Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "Centrality estimator (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3)"};
  Configurable<bool> cfgPerformPidSelectionInIts{"cfgPerformPidSelectionInIts", false, "Perform PID selections in ITS"};
  Configurable<bool> cfgPerformPidSelectionInTpc{"cfgPerformPidSelectionInTpc", false, "Perform PID selections in TPC"};
  Configurable<bool> cfgPerformPidSelectionInTof{"cfgPerformPidSelectionInTof", false, "Perform PID selections in TOF"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {nuclei::betheBlochDefault[0], nuclei::Species::kNspecies, 6, nuclei::names, nuclei::betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<int>> cfgUseCentralTpcCalibration{"cfgUseCentralTpcCalibration", {nuclei::useCentralTpcCalibrationDefault[0], nuclei::Species::kNspecies, 1, nuclei::names, {"UseCentralTpcCalibration"}}, "Use central TPC calibration"};
  Configurable<LabeledArray<double>> cfgDownscalingFactor{"cfgDownscalingFactor", {nuclei::DownscalingDefault[0], nuclei::Species::kNspecies, 1, nuclei::names, {"DownscalingFactor"}}, "Save to the AO2D with a downscaling factor"};

  Configurable<LabeledArray<int>> cfgUseTrackTuner{"cfgUseTrackTuner", {nuclei::useTrackTuner[0], nuclei::Species::kNspecies, 1, nuclei::names, {"UseTrckTuner"}}, "Use Track Tuner"};
  Configurable<std::string> cfgTrackTunerParams{"cfgTrackTunerParams", "debugInfo=0|updateTrackDCAs=1|updateTrackCovMat=1|updateCurvature=0|updateCurvatureIU=0|updatePulls=0|isInputFileFromCCDB=1|pathInputFile=Users/m/mfaggin/test/inputsTrackTuner/pp2023/pass4/correct_names|nameInputFile=trackTuner_DataLHC23hPass4_McLHC23k4g.root|pathFileQoverPt=Users/h/hsharma/qOverPtGraphs|nameFileQoverPt=D0sigma_Data_removal_itstps_MC_LHC22b1b.root|usePvRefitCorrections=0|qOverPtMC=-1.|qOverPtData=-1.|nPhiBins=1|autoDetectDcaCalib=false", "TrackTuner parameter initialization (format: <name>=<value>|<name>=<value>)"};
  Configurable<int> cfgTrackTunerConfigSource{"cfgTrackTunerConfigSource", aod::track_tuner::InputString, "1: input string; 2: TrackTuner Configurables"};
  ConfigurableAxis cfgAxisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  Configurable<float> cfgRapidityMin{"cfgRapidityMin", -1., "Minimum rapidity value"};
  Configurable<float> cfgRapidityMax{"cfgRapidityMax", 1., "Maximum rapidity value"};
  Configurable<float> cfgRapidityCenterMass{"cfgRapidityCenterMass", 0.0f, "Center of mass rapidity"};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutTpcMom{"cfgCutTpcMom", 0.2f, "Minimum TPC momentum for tracks"};
  Configurable<float> cfgCutNclusITS{"cfgCutNclusITS", 5, "Minimum number of ITS clusters"};
  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};
  Configurable<float> cfgCutNclusCrossedRowsTPC{"cfgCutNclusCrossedRowsTPC", 70, "Minimum number of TPC clusters crossed rows"};
  Configurable<float> cfgCutChi2PerClusterTPC{"cfgCutChi2PerClusterTPC", 4.f, "Maximum chi2 per TPC cluster"};
  Configurable<float> cfgCutChi2PerClusterITS{"cfgCutChi2PerClusterITS", 36.f, "Maximum chi2 per ITS cluster"};

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
      {"hFailCentrality", "0: all the times the centrality filling function is called - 1: each time it fails ; Bool", {HistType::kTH1F, {{2, -0.5, 1.50}}}},
      {"hTrackTunedTracks", "", {HistType::kTH1F, {{1, 0.5, 1.5}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};
  std::shared_ptr<TH1> mHistFailCentrality = mHistograms.get<TH1>(HIST("hFailCentrality"));
  std::shared_ptr<TH1> mHistTrackTunedTracks = mHistograms.get<TH1>(HIST("hTrackTunedTracks"));

  std::vector<int> mSpeciesToProcess;
  Produces<aod::NucleiTableRed> mNucleiTableRed;
  Produces<aod::NucleiTableExt> mNucleiTableExt;

  std::vector<nuclei::SlimCandidate> mNucleiCandidates;
  std::vector<int> mFilledMcParticleIds;

  TrackTuner mTrackTuner;
  o2::base::Propagator::MatCorrType mMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  std::array<float, 2> mDcaInfo;
  o2::dataformats::DCA mDcaInfoCov;
  o2::dataformats::VertexBase mVtx;
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
      if (cfgSpeciesToProcess->get(iSpecies) == 1)
        mSpeciesToProcess.emplace_back(iSpecies);
    }

    static_for<0, nuclei::kNspecies - 1>([&](auto iSpecies) {
      constexpr int kSpeciesCt = decltype(iSpecies)::value;
      const int kSpeciesRt = kSpeciesCt;

      if (std::find(mSpeciesToProcess.begin(), mSpeciesToProcess.end(), kSpeciesCt) == mSpeciesToProcess.end())
        return;

      float tpcBetheBlochParams[6];
      for (int iParam = 0; iParam < 6; iParam++) {
        tpcBetheBlochParams[iParam] = cfgBetheBlochParams->get(kSpeciesRt, iParam);
      }

      nuclei::createHistogramRegistryNucleus<kSpeciesCt>(mHistograms);
      mHistograms.add(fmt::format("{}/hTrackQuality", nuclei::cNames[kSpeciesRt]).c_str(), (fmt::format("{} track quality;", nuclei::cNames[kSpeciesRt]) + std::string("#it{p}_{T} / #it{Z} (GeV/#it{c}); Selection step; Counts")).c_str(), o2::framework::HistType::kTH2D, {{400, -10.0f, 10.0f}, {trackQuality::kNtrackQuality, -0.5f, static_cast<float>(trackQuality::kNtrackQuality) - 0.5f}});
      for (size_t iSel = 0; iSel < trackQuality::kNtrackQuality; iSel++) {
        mHistograms.get<TH2>(HIST(nuclei::cNames[kSpeciesRt]) + HIST("/hTrackQuality"))->GetYaxis()->SetBinLabel(iSel + 1, trackQualityLabels[iSel].c_str());
      }

      if (cfgUseCentralTpcCalibration->get(static_cast<uint32_t>(kSpeciesRt), static_cast<uint32_t>(0)) == 0) {
        mPidManagers[kSpeciesRt] = nuclei::PidManager(kSpeciesRt, tpcBetheBlochParams);
      } else {
        mPidManagers[kSpeciesRt] = nuclei::PidManager(kSpeciesRt);
      }
    });

    /// TrackTuner initialization
    bool anyTrackTuner = false;
    for (int iSpecies = 0; iSpecies < static_cast<int>(nuclei::Species::kNspecies); iSpecies++) {
      anyTrackTuner = anyTrackTuner || cfgUseTrackTuner->get(iSpecies);
    }

    if (anyTrackTuner) {
      std::string outputStringParams = "";
      switch (cfgTrackTunerConfigSource) {
        case aod::track_tuner::InputString:
          outputStringParams = mTrackTuner.configParams(cfgTrackTunerParams);
          break;
        case aod::track_tuner::Configurables:
          outputStringParams = mTrackTuner.configParams();
          break;

        default:
          LOG(fatal) << "TrackTuner configuration source not defined. Fix it! (Supported options: input string (1); Configurables (2))";
          break;
      }

      if (!mTrackTuner.autoDetectDcaCalib) {
        mTrackTuner.getDcaGraphs();
        mHistTrackTunedTracks->SetTitle(outputStringParams.c_str());
      }
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
    o2::base::Propagator::Instance()->setMatLUT(nuclei::lut);
    mBz = static_cast<float>(grpmag->getNominalL3Field());
    LOGF(info, "Retrieved GRP for timestamp %ull (%i) with magnetic field of %1.2f kZG", timestamp, mRunNumber, mBz);
  }

  template <int iSpecies, typename Ttrack>
  bool trackSelection(const Ttrack& track)
  {
    mHistograms.fill(HIST(nuclei::cNames[iSpecies]) + HIST("/hTrackQuality"), track.sign() * track.pt(), trackQuality::kNoCuts);

    if (std::abs(track.eta()) > cfgCutEta)
      return false;
    mHistograms.fill(HIST(nuclei::cNames[iSpecies]) + HIST("/hTrackQuality"), track.sign() * track.pt(), trackQuality::kEtaCut);

    if (track.tpcInnerParam() < cfgCutTpcMom)
      return false;
    mHistograms.fill(HIST(nuclei::cNames[iSpecies]) + HIST("/hTrackQuality"), track.sign() * track.pt(), trackQuality::kPtpcCut);

    if (track.itsNCls() < cfgCutNclusITS)
      return false;
    mHistograms.fill(HIST(nuclei::cNames[iSpecies]) + HIST("/hTrackQuality"), track.sign() * track.pt(), trackQuality::kNclsItsCut);

    if (track.tpcNClsFound() < cfgCutNclusTPC)
      return false;
    mHistograms.fill(HIST(nuclei::cNames[iSpecies]) + HIST("/hTrackQuality"), track.sign() * track.pt(), trackQuality::kNclsTpcCut);

    if (track.tpcNClsCrossedRows() < cfgCutNclusCrossedRowsTPC)
      return false;
    mHistograms.fill(HIST(nuclei::cNames[iSpecies]) + HIST("/hTrackQuality"), track.sign() * track.pt(), trackQuality::kNCrossedRowsCut);

    if (track.tpcChi2NCl() > cfgCutChi2PerClusterTPC)
      return false;
    mHistograms.fill(HIST(nuclei::cNames[iSpecies]) + HIST("/hTrackQuality"), track.sign() * track.pt(), trackQuality::kTpcChi2Cut);

    if (track.itsChi2NCl() > cfgCutChi2PerClusterITS)
      return false;
    mHistograms.fill(HIST(nuclei::cNames[iSpecies]) + HIST("/hTrackQuality"), track.sign() * track.pt(), trackQuality::kItsChi2Cut);

    return true;
  }

  template <int iSpecies, typename Ttrack, typename Tcollision>
  bool pidSelection(const Ttrack& track, const Tcollision& collision)
  {
    constexpr int kIndex = iSpecies;
    if (!nuclei::checkSpeciesValidity(kIndex))
      std::runtime_error("species contains invalid nucleus kIndex");

    float centrality = nuclei::getCentrality(collision, cfgCentralityEstimator);
    float nsigmaTPC = mPidManagers[kIndex].getNSigmaTPC(track);
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaTPC_preselectionVsCentrality"), track.pt() * track.sign(), nsigmaTPC, centrality);
    if (std::abs(nsigmaTPC) > cfgNsigmaTPC->get(kIndex, 1) && cfgPerformPidSelectionInTpc)
      return false;
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaTPCVsCentrality"), track.pt() * track.sign(), nsigmaTPC, centrality);

    float nsigmaITS = mPidManagers[kIndex].getNSigmaITS(track);
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaITS_preselectionVsCentrality"), track.sign() * track.pt(), nsigmaITS, centrality);
    // add nsigmaITS cut ?
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaITSVsCentrality"), track.sign() * track.pt(), nsigmaITS, centrality);

    float nsigmaTOF = mPidManagers[kIndex].getNSigmaTOF(track);
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaTOF_preselectionVsCentrality"), track.sign() * track.pt(), nsigmaTOF, centrality);
    if (std::abs(nsigmaTOF) > cfgNsigmaTOF->get(kIndex, 1) && track.hasTOF() && cfgPerformPidSelectionInTof)
      return false;
    mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3NsigmaTOFVsCentrality"), track.sign() * track.pt(), nsigmaTOF, centrality);

    return true;
  }

  template <typename Tparticle>
  void fillNucleusFlagsPdgsMc(const Tparticle& particle, nuclei::SlimCandidate& candidate)
  {
    candidate.pdgCode = particle.pdgCode();
    candidate.mcProcess = particle.getProcess();

    if (particle.has_mothers()) {
      for (const auto& motherparticle : particle.template mothers_as<aod::McParticles>()) {
        candidate.motherPdgCode = motherparticle.pdgCode();
      }
    }

    if (particle.isPhysicalPrimary()) {
      candidate.flags |= nuclei::Flags::kIsPhysicalPrimary;

      ///<  heavy flavour mother
      /*if (particle.has_mothers()) {
        for (const auto& motherparticle : particle.mothers_as<aod::McParticles>()) {
          if (std::find(nuclei::hfMothCodes.begin(), nuclei::hfMothCodes.end(), std::abs(motherparticle.pdgCode())) != nuclei::hfMothCodes.end()) {
            flags |= kIsSecondaryFromWeakDecay;
            motherPdgCode = motherparticle.pdgCode();
            break;
          }
        }
      }*/

    } else if (particle.getProcess() == TMCProcess::kPDecay) {
      ///<  assuming that strong decays are included in the previous step
      candidate.flags |= nuclei::Flags::kIsSecondaryFromWeakDecay;
    } else {
      candidate.flags |= nuclei::Flags::kIsSecondaryFromMaterial;
    }
  }

  template <typename Tcollision, typename Ttrack>
  void fillNucleusFlagsPdgs(const int iSpecies, const Tcollision& collision, const Ttrack& track, nuclei::SlimCandidate& candidate)
  {
    candidate.flags = static_cast<uint16_t>((track.pidForTracking() & 0xF) << 12);

    switch (iSpecies) {
      case nuclei::Species::kPr:
        candidate.flags |= nuclei::Flags::kProton;
        break;
      case nuclei::Species::kDe:
        candidate.flags |= nuclei::Flags::kDeuteron;
        break;
      case nuclei::Species::kTr:
        candidate.flags |= nuclei::Flags::kTriton;
        break;
      case nuclei::Species::kHe:
        candidate.flags |= nuclei::Flags::kHe3;
        break;
      case nuclei::Species::kAl:
        candidate.flags |= nuclei::Flags::kHe4;
        break;
      default:
        candidate.flags |= 0;
        break;
    }

    if (track.hasTOF())
      candidate.flags |= nuclei::Flags::kHasTOF;

    if (track.hasTRD())
      candidate.flags |= nuclei::Flags::kHasTRD;

    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      candidate.flags |= nuclei::Flags::kITSrof;
  }

  template <typename Tparticle>
  void fillNucleusGeneratedVariables(const Tparticle& particle, nuclei::SlimCandidate& candidate)
  {
    candidate.ptGenerated = particle.pt() * (particle.pdgCode() > 0 ? 1.f : -1.f);
    candidate.etaGenerated = particle.eta();
    candidate.yGenerated = particle.y();
    candidate.phiGenerated = particle.phi();
  }

  template <const bool isMc, typename Tcollision, typename Ttrack>
  void fillDcaInformation(const int iSpecies, const Tcollision& collision, const Ttrack& track, nuclei::SlimCandidate& candidate, const aod::McParticles::iterator& particle)
  {

    const o2::math_utils::Point3D<float> collisionVertex{collision.posX(), collision.posY(), collision.posZ()};

    mDcaInfoCov.set(999, 999, 999, 999, 999);
    setTrackParCov(track, mTrackParCov);
    mTrackParCov.setPID(track.pidForTracking());

    if constexpr (isMc) {
      if (track.has_mcParticle() && cfgUseTrackTuner->get(iSpecies)) {
        mHistTrackTunedTracks->Fill(1.);
        mTrackTuner.tuneTrackParams(particle, mTrackParCov, mMatCorr, &mDcaInfoCov, mHistTrackTunedTracks);
      }
    } else {
      mMatCorr = static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value);
    }

    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCov, 2.f, mMatCorr, &mDcaInfoCov);

    candidate.DCAxy = mDcaInfoCov.getY();
    candidate.DCAz = mDcaInfoCov.getZ();
  }

  template <const bool isMc, typename Tcollision, typename Ttrack>
  nuclei::SlimCandidate fillCandidate(const int iSpecies, Tcollision const& collision, Ttrack const& track)
  {
    if (!nuclei::checkSpeciesValidity(iSpecies))
      std::runtime_error("species contains invalid nucleus index");

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
                                       .yGenerated = 0.f,
                                       .phiGenerated = 0.f,
                                       .centrality = nuclei::getCentrality(collision, cfgCentralityEstimator, mHistFailCentrality),
                                       .mcProcess = TMCProcess::kPNoProcess,
                                       .nsigmaTpc = mPidManagers[iSpecies].getNSigmaTPC(track),
                                       .nsigmaTof = mPidManagers[iSpecies].getNSigmaTOF(track)};

    fillNucleusFlagsPdgs(iSpecies, collision, track, candidate);

    aod::McParticles::iterator particle;

    if constexpr (isMc) {
      if (track.has_mcParticle()) {

        particle = track.mcParticle();
        fillNucleusFlagsPdgsMc(particle, candidate);
        fillNucleusGeneratedVariables(particle, candidate);
      }
    }

    fillDcaInformation<isMc>(iSpecies, collision, track, candidate, particle);

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
    if (!nuclei::checkSpeciesValidity(kIndex))
      std::runtime_error("species contains invalid nucleus kIndex");

    if (isGenerated) {
      const float ptGenerated = (kIndex == nuclei::Species::kPr || kIndex == nuclei::Species::kDe || kIndex == nuclei::Species::kTr) ? candidate.ptGenerated : candidate.ptGenerated / 2.f;
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/hPtGenerated"), ptGenerated);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h2PtVsCentralityGenerated"), ptGenerated, candidate.centrality);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3PtVsRapidityVsCentralityGenerated"), ptGenerated, candidate.yGenerated, candidate.centrality);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h3PhiVsEtaVsCentralityGenerated"), candidate.phiGenerated, candidate.etaGenerated, candidate.centrality);
    } else {
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/hPtReconstructed"), candidate.pt);
      mHistograms.fill(HIST(nuclei::cNames[kIndex]) + HIST("/h2PtVsCentralityReconstructed"), candidate.pt, candidate.centrality);
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

    gRandom->SetSeed(67);
    mNucleiCandidates.clear();
    mFilledMcParticleIds.clear();

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!nuclei::eventSelection(collision, mHistograms, cfgEventSelections, cfgCutVertex))
      return;

    bool anyTrackTuner = false;
    for (int iSpecies = 0; iSpecies < static_cast<int>(nuclei::Species::kNspecies); iSpecies++) {
      anyTrackTuner = anyTrackTuner || cfgUseTrackTuner->get(iSpecies);
    }
    if (anyTrackTuner && mTrackTuner.autoDetectDcaCalib && !mTrackTuner.areGraphsConfigured) {

      mTrackTuner.setRunNumber(mRunNumber);

      /// setup the "auto-detected" path based on the run number
      mTrackTuner.getPathInputFileAutomaticFromCCDB();
      mHistTrackTunedTracks->SetTitle(mTrackTuner.outputString.c_str());
      mTrackTuner.getDcaGraphs();
    }

    auto tracksThisCollision = tracks.sliceBy(mTracksPerCollision, collision.globalIndex());
    tracksThisCollision.bindExternalIndices(&tracks);

    for (const auto& track : tracks) {

      static_for<0, nuclei::kNspecies - 1>([&](auto iSpecies) {
        constexpr int kSpeciesCt = decltype(iSpecies)::value;
        const int kSpeciesRt = kSpeciesCt;

        if (std::find(mSpeciesToProcess.begin(), mSpeciesToProcess.end(), kSpeciesRt) == mSpeciesToProcess.end())
          return;

        if (!track.has_mcParticle())
          return;

        const auto& particle = track.mcParticle();
        if (cfgDoCheckPdgCode) {
          if (std::abs(particle.pdgCode()) != nuclei::pdgCodes[kSpeciesRt])
            return;
        }

        if (cfgDownscalingFactor->get(kSpeciesRt) < 1.) {
          if ((gRandom->Uniform()) > cfgDownscalingFactor->get(kSpeciesRt))
            return;
        }

        if ((particle.y() - cfgRapidityCenterMass) < cfgRapidityMin || (particle.y() - cfgRapidityCenterMass) > cfgRapidityMax)
          return;

        if (cfgFillOnlyPhysicalPrimaries && !particle.isPhysicalPrimary())
          return;

        mHistograms.fill(HIST(nuclei::cNames[kSpeciesCt]) + HIST("/hTrackSelections"), nuclei::trackSelection::kNoCuts);
        if (!trackSelection<kSpeciesRt>(track))
          return;
        mHistograms.fill(HIST(nuclei::cNames[kSpeciesCt]) + HIST("/hTrackSelections"), nuclei::trackSelection::kTrackCuts);

        if (!pidSelection<kSpeciesRt>(track, collision))
          return;
        mHistograms.fill(HIST(nuclei::cNames[kSpeciesCt]) + HIST("/hTrackSelections"), nuclei::trackSelection::kPidCuts);

        nuclei::SlimCandidate candidate;
        candidate = fillCandidate</*isMc*/ true>(kSpeciesCt, collision, track);

        mNucleiCandidates.emplace_back(candidate);
        mFilledMcParticleIds.emplace_back(particle.globalIndex());
        dispatchFillHistograms</*isGenerated*/ true>(kSpeciesRt, candidate);
        dispatchFillHistograms</*isGenerated*/ false>(kSpeciesRt, candidate);
      });
    }

    const int mcCollisionId = collision.mcCollisionId();
    auto mcParticlesThisCollision = mcParticles.sliceBy(mMcParticlesPerCollision, mcCollisionId);
    mcParticlesThisCollision.bindExternalIndices(&mcParticles);

    for (const auto& particle : mcParticlesThisCollision) {

      if (std::find(mFilledMcParticleIds.begin(), mFilledMcParticleIds.end(), particle.globalIndex()) != mFilledMcParticleIds.end())
        continue;

      if (cfgFillOnlyPhysicalPrimaries && !particle.isPhysicalPrimary())
        continue;

      if ((particle.y() - cfgRapidityCenterMass) < cfgRapidityMin || (particle.y() - cfgRapidityCenterMass) > cfgRapidityMax)
        continue;

      int iSpecies = nuclei::getSpeciesFromPdg(particle.pdgCode());
      if (std::find(mSpeciesToProcess.begin(), mSpeciesToProcess.end(), iSpecies) == mSpeciesToProcess.end())
        continue;

      if (cfgDownscalingFactor->get(iSpecies) < 1.) {
        if ((gRandom->Uniform()) > cfgDownscalingFactor->get(iSpecies))
          return;
      }

      nuclei::SlimCandidate candidate;
      candidate.centrality = nuclei::getCentrality(collision, cfgCentralityEstimator, mHistFailCentrality);
      fillNucleusFlagsPdgsMc(particle, candidate);
      fillNucleusGeneratedVariables(particle, candidate);

      mNucleiCandidates.emplace_back(candidate);
      mFilledMcParticleIds.emplace_back(particle.globalIndex());
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
        candidate.ptGenerated,
        candidate.mcProcess,
        candidate.pdgCode,
        candidate.motherPdgCode);
      mNucleiTableExt(
        candidate.nsigmaTpc,
        candidate.nsigmaTof);
    }
  }
  PROCESS_SWITCH(nucleiQC, processMc, "Mc analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<nucleiQC>(cfgc)};
}
