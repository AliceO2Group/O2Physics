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

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <array>

#include "Math/Vector4D.h"

#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/Tools/TrackTuner.h"
#include "Common/Core/RecoDecay.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/Track.h"

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFSlimNucleiTables.h"

#include "TRandom3.h"

#include "PWGLF/Utils/nucleiUtils.h"

using namespace o2;
using namespace o2::framework;


struct nucleiQC {

//using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksDCA, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCDe, aod::pidTOFDe, aod::pidTPCTr, aod::pidTOFTr, aod::pidTPCHe, aod::pidTOFHe, aod::pidTPCAl, aod::pidTOFAl>;
  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksDCA, aod::TracksExtra, aod::pidTPCPr, aod::pidTPCDe, aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl>;
//using TrackCandidatesMC = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksDCA, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCDe, aod::pidTOFDe, aod::pidTPCTr, aod::pidTOFTr, aod::pidTPCHe, aod::pidTOFHe, aod::pidTPCAl, aod::pidTOFAl, aod::McTrackLabels>;
  using TrackCandidatesMC = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksDCA, aod::TracksExtra, aod::pidTPCPr, aod::pidTPCDe, aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl, aod::McTrackLabels>;
  using Collision = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>::iterator;

  Configurable<bool> cfgFillTable{"cfgFillTable", true, "Fill output tree"};
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

  HistogramRegistry mQaHistograms {
    "QA",
    {
      {"hEventSelections", "Event selections; Selection step; Counts", {HistType::kTH1D, {{nuclei::evSel::kNevSels + 1, -0.5f, static_cast<float>(nuclei::evSel::kNevSels) + 0.5f}}}},
      {"hVtxZBefore", "Vertex distribution in Z before selections;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
      {"hVtxZ", "Vertex distribution in Z;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true
  };
  std::array<HistogramRegistry, static_cast<int>(nuclei::Species::kNspecies)> mNucleiHistogramRegistries;
  std::vector<int> mSpeciesToProcess;
  Produces<aod::NucleiTableRed> mNucleiTableRed;
  
  std::vector<nuclei::SlimCandidate> mNucleiCandidates;
  std::vector<int> mFilledMcParticleIds;
  
  std::array<nuclei::PidManager, static_cast<int>(nuclei::Species::kNspecies)> mPidManagers;

  void init(o2::framework::InitContext&) {

    LOG(info) << "Check init";

    for (int iSel = 0; iSel < nuclei::evSel::kNevSels; iSel++) {
      mQaHistograms.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(iSel + 1, nuclei::eventSelectionLabels[iSel].c_str());
    }

    LOG(info) << " check after hist ev sel";

    for (int iSpecies = 0; iSpecies < static_cast<int>(nuclei::Species::kNspecies); iSpecies++) {
      if (cfgSpeciesToProcess->get(iSpecies) == 1) {
        mSpeciesToProcess.emplace_back(iSpecies);
      }
    }

    LOG(info) << " check after species to be processed were saved";

    for (auto iSpecies: mSpeciesToProcess) {

      float tpcBetheBlochParams[6];
      for (int iParam = 0; iParam < 6; iParam++) {
        tpcBetheBlochParams[iParam] = cfgBetheBlochParams->get(iSpecies, iParam);
      }

      LOG(info) << " check after bb params, ispecies: " << iSpecies;

      mNucleiHistogramRegistries[iSpecies] = nuclei::createHistogramRegistryNucleus(iSpecies);
      
      if (cfgUseCentralTpcCalibration->get(uint32_t(iSpecies), uint32_t(0)) == 0) {
        mPidManagers[iSpecies] = nuclei::PidManager(nuclei::Species(iSpecies), tpcBetheBlochParams);
      } else {
        mPidManagers[iSpecies] = nuclei::PidManager(nuclei::Species(iSpecies));
      }

      LOG(info) << " check after pidManagers, ispecies: " << iSpecies;
    }

    LOG(info) << " check end init";
  }

  template <typename Ttrack>
  bool trackSelection(const Ttrack& track) {
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
  
  template <typename Ttrack, typename Tcollision>
  bool pidSelection(const int iSpecies, const Ttrack& track, const Tcollision& collision) 
  { 
    if (!nuclei::checkSpeciesValidity(iSpecies)) {
      std::runtime_error("species contains invalid nucleus index");
    }

    const float centrality = nuclei::getCentrality(collision, cfgCentralityEstimator);
    const float nsigmaTPC = mPidManagers[iSpecies].getNSigmaTPC(track);
    mNucleiHistogramRegistries[iSpecies].fill(HIST("h2NsigmaTPC_preselectionVsCentrality"), track.pt(), nsigmaTPC, centrality);
    if (std::abs(nsigmaTPC) > cfgNsigmaTPC->get(iSpecies, 1)) {
      return false;
    }
    mNucleiHistogramRegistries[iSpecies].fill(HIST("h2NsigmaTPCVsCentrality"), track.pt(), nsigmaTPC, centrality);

    const float nsigmaITS = mPidManagers[iSpecies].getNSigmaITS(track);
    mNucleiHistogramRegistries[iSpecies].fill(HIST("h2NSigmaITS_preselectionVsCentrality"), track.sign() * track.pt(), nsigmaITS, centrality);
    // add nsigmaITS cut ?
    mNucleiHistogramRegistries[iSpecies].fill(HIST("h2NSigmaITSVsCentrality"), track.sign() * track.pt(), nsigmaITS, centrality);

    //const float nsigmaTOF = mPidManagers[iSpecies].getNSigmaTOF(track);
    //mNucleiHistogramRegistries[iSpecies].fill(HIST("h2NSigmaTOF_preselectionVsCentrality"), track.pt(), nsigmaTOF, centrality);
    //if (std::abs(nsigmaTOF) > cfgNsigmaTOF->get(iSpecies, 1)) {
    //  return false;
    //}
    //mNucleiHistogramRegistries[iSpecies].fill(HIST("h2NSigmaTOF"), track.pt(), nsigmaTOF);

    return true;
  }

  template <typename Tparticle>
  void fillNucleusFlagsPdgsMc(const Tparticle& particle, nuclei::SlimCandidate& candidate) 
  {
    candidate.pdgCode = particle.pdgCode();

    if (particle.isPhysicalPrimary()) {
        candidate.flags |= nuclei::Flags::kIsPhysicalPrimary;
        
        // heavy flavour mother 
        //if (particle.has_mothers()) {
        //  for (auto& motherparticle : particle.mothers_as<aod::McParticles>()) {
        //    if (std::find(nuclei::hfMothCodes.begin(), nuclei::hfMothCodes.end(), std::abs(motherparticle.pdgCode())) != nuclei::hfMothCodes.end()) {
        //      flags |= kIsSecondaryFromWeakDecay;
        //      motherPdgCode = motherparticle.pdgCode();
        //      break;
        //    }
        //  }
        //}

      } else if (particle.has_mothers()) {
        candidate.flags |= nuclei::Flags::kIsSecondaryFromWeakDecay;
        for (auto& motherparticle : particle.template mothers_as<aod::McParticles>()) {
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
    candidate.flags |= iSpecies == nuclei::Species::kPr ? nuclei::Flags::kProton :
                       iSpecies == nuclei::Species::kDe ? nuclei::Flags::kDeuteron :
                       iSpecies == nuclei::Species::kTr ? nuclei::Flags::kTriton :
                       iSpecies == nuclei::Species::kHe ? nuclei::Flags::kHe3 :
                       iSpecies == nuclei::Species::kAl ? nuclei::Flags::kHe4 : 0;
    
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
                                       //.beta = mPidManagers[iSpecies].getBetaTOF(track),
                                       .beta = 0.f,
                                       .DCAxy = track.dcaXY(),
                                       .DCAz = track.dcaZ(),
                                       .flags = 0,
                                       .pdgCode = 0,
                                       .motherPdgCode = 0,
                                       .ptGenerated = 0.f, // to be filled for mc
                                       .etaGenerated = 0.f,
                                       .phiGenerated = 0.f,
                                       .centrality = nuclei::getCentrality(collision, cfgCentralityEstimator)
                                       };
    
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

  template <const bool isGenerated>
  void fillHistograms(const int iSpecies, const nuclei::SlimCandidate& candidate)
  {
    if (!nuclei::checkSpeciesValidity(iSpecies)) {
      std::runtime_error("species contains invalid nucleus index");
    }

    if (isGenerated) {
      const float ptGenerated = (iSpecies == nuclei::Species::kPr || iSpecies == nuclei::Species::kDe || iSpecies == nuclei::Species::kTr) ? candidate.ptGenerated : candidate.ptGenerated / 2.f;
      mNucleiHistogramRegistries[iSpecies].fill(HIST("hPtGenerated"), ptGenerated);
      mNucleiHistogramRegistries[iSpecies].fill(HIST("h3PtVsEtaVsCentralityGenerated"), ptGenerated, candidate.etaGenerated, candidate.centrality);
      mNucleiHistogramRegistries[iSpecies].fill(HIST("hPhiPhiVsEtaVsCentralityGenerated"), candidate.phiGenerated, candidate.etaGenerated, candidate.centrality);
    } else {
      mNucleiHistogramRegistries[iSpecies].fill(HIST("hPt"), candidate.pt);
      mNucleiHistogramRegistries[iSpecies].fill(HIST("h3PtVsEtaVsCentralityReconstructed"), candidate.pt, candidate.eta, candidate.centrality);
      mNucleiHistogramRegistries[iSpecies].fill(HIST("hPhiPhiVsEtaVsCentralityReconstructed"), candidate.phi, candidate.eta, candidate.centrality);
      mNucleiHistogramRegistries[iSpecies].fill(HIST("h3DCAxyVsPtVsCentrality"), candidate.DCAxy, candidate.pt, candidate.centrality);
      mNucleiHistogramRegistries[iSpecies].fill(HIST("h3DCAzVsPtVsCentrality"), candidate.DCAz, candidate.pt, candidate.centrality);
      mNucleiHistogramRegistries[iSpecies].fill(HIST("h3BetaVsPtVsCentrality"), candidate.beta, candidate.pt, candidate.centrality);
      mNucleiHistogramRegistries[iSpecies].fill(HIST("h3dEdxVsPVsCentrality"), candidate.TPCsignal, candidate.pt, candidate.centrality);
      mNucleiHistogramRegistries[iSpecies].fill(HIST("h3ClusterSizeVsPtVsCentrality"), mPidManagers[iSpecies].getClusterSizeCosLambdaITS(candidate.clusterSizesITS, candidate.eta), candidate.pt, candidate.centrality);
    }
  }

  void processMc(const Collision& collision, const TrackCandidatesMC& tracks, const aod::BCsWithTimestamps&, const aod::McParticles& mcParticles)
  {
    mNucleiCandidates.clear();
    mFilledMcParticleIds.clear();

    LOG(info) << "Running";

    for (const auto iSpecies: mSpeciesToProcess) {
      mPidManagers[iSpecies].initMcResponseITS();
    } 

    if (!nuclei::eventSelection(collision, mQaHistograms, cfgEventSelections, cfgCutVertex)) {
      return;
    }

    for (const auto& track: tracks)
    {
      for (const auto iSpecies: mSpeciesToProcess) {
        
        mNucleiHistogramRegistries[iSpecies].fill(HIST("hTrackSelections"), nuclei::trackSelection::kNoCuts);
        if (!trackSelection(track)) {
          continue;
        }

        mNucleiHistogramRegistries[iSpecies].fill(HIST("hTrackSelections"), nuclei::trackSelection::kTrackCuts);
        if (!pidSelection(iSpecies, track, collision)) {
          continue;
        }

        mNucleiHistogramRegistries[iSpecies].fill(HIST("hTrackSelections"), nuclei::trackSelection::kPidCuts);
        
        nuclei::SlimCandidate candidate;
        if (track.has_mcParticle()) {
          const auto& particle = track.mcParticle();
          if (particle.y() - cfgRapidityCenterMass < cfgRapidityMin || particle.y() - cfgRapidityCenterMass > cfgRapidityMax) {
            continue;
          }
          candidate = fillCandidate</*isMc*/true>(iSpecies, collision, track);
          fillHistograms</*isGenerated*/true>(iSpecies, candidate);
        } else {
          candidate = fillCandidate</*isMc*/true>(iSpecies, collision, track);
        }
        
        fillHistograms</*isGenerated*/false>(iSpecies, candidate);
      }
    }

    for (const auto& particle: mcParticles)
    {
      if (std::find(mFilledMcParticleIds.begin(), mFilledMcParticleIds.end(), particle.globalIndex()) != mFilledMcParticleIds.end()) {
        continue;
      }

      nuclei::SlimCandidate candidate;
      fillNucleusFlagsPdgsMc(particle, candidate);
      fillNucleusGeneratedVariables(particle, candidate);
      mNucleiCandidates.emplace_back(candidate);

      const int iSpecies = nuclei::getSpeciesFromPdg(particle.pdgCode());
      fillHistograms</*isGenerated*/true>(iSpecies, candidate);
    }

    if (!cfgFillTable) {
      return;
    }

    for (const auto& candidate: mNucleiCandidates)
    {
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
        candidate.motherPdgCode
      );
    }
  }
  PROCESS_SWITCH(nucleiQC, processMc, "Mc analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<nucleiQC>(cfgc)};
}
