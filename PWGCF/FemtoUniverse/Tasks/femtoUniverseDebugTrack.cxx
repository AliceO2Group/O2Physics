// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoUniverseDebugTrack.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for tracks
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int NCuts = 5;
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[1][NCuts] = {{4.05f, 0.75f, 3.f, 3.f, 100.f}};

} // namespace

struct FemtoUniverseDebugTrack {
  SliceCache cache;

  Configurable<LabeledArray<float>> confCutTable{"confCutTable", {cutsTable[0], NCuts, cutNames}, "Particle selections"};
  Configurable<int> confNspecies{"confNspecies", 2, "Number of particle spieces with PID info"};

  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confPDGCodePartOne{"confPDGCodePartOne", 2212, "Particle 1 - PDG code"};
    Configurable<bool> confIsTrackIdentified{"confIsTrackIdentified", true, "Enable PID for the track"};
    Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<int> confTrackChoicePartOne{"confTrackChoicePartOne", 1, "Type of particle (track1): {0:Proton, 1:Pion, 2:Kaon}"};
    Configurable<int> confTrackSign{"confTrackSign", 1, "Track sign"};
  } trackonefilter;

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confNsigmaCombinedKaon{"confNsigmaCombinedKaon", 3.0, "TPC and TOF Kaon Sigma (combined) for momentum > 0.4"};
    Configurable<float> confNsigmaTPCKaon{"confNsigmaTPCKaon", 3.0, "TPC Kaon Sigma for momentum < 0.4"};
    Configurable<float> confNsigmaCombinedProton{"confNsigmaCombinedProton", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > 0.5"};
    Configurable<float> confNsigmaTPCProton{"confNsigmaTPCProton", 3.0, "TPC Proton Sigma for momentum < 0.5"};
    Configurable<float> confNsigmaCombinedPion{"confNsigmaCombinedPion", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > 0.5"};
    Configurable<float> confNsigmaTPCPion{"confNsigmaTPCPion", 3.0, "TPC Pion Sigma for momentum < 0.5"};
  } generalPIDcuts;

  // Configurable<uint32_t> confCutPartOne{"confCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
  Configurable<int> confPIDPartOne{"confPIDPartOne", 1, "Particle 1 - Read from cutCulator"};
  Configurable<std::vector<float>> confTrkPIDnSigmaMax{"confTrkPIDnSigmaMax", std::vector<float>{3.5f, 3.f, 2.5f}, "This configurable needs to be the same as the one used in the producer task"};
  ConfigurableAxis confTempFitVarBins{"confDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarpTBins{"confTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  FemtoUniverseTrackSelection trackCuts;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Partition<FemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == int8_t(trackonefilter.confTrackSign)); // && ((aod::femtouniverseparticle::cut & confCutPartOne) == confCutPartOne);
  Preslice<FemtoFullParticles> perColReco = aod::femtouniverseparticle::fdCollisionId;

  using FemtoFullParticlesMC = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Partition<FemtoFullParticlesMC> partsOneMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == int8_t(trackonefilter.confTrackSign)); // && ((aod::femtouniverseparticle::cut & confCutPartOne) == confCutPartOne);
  Preslice<FemtoFullParticlesMC> perColGen = aod::femtouniverseparticle::fdCollisionId;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack> trackHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPartOne;
  std::vector<float> kNsigma;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  // PID for protons
  bool IsProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // confNsigmaTPCProton -> TPC Kaon Sigma for momentum < 0.5
    // confNsigmaCombinedProton -> TPC and TOF Kaon Sigma (combined) for momentum > 0.5

    if (mom < 0.5) {
      if (std::abs(nsigmaTPCPr) < generalPIDcuts.confNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 0.5) {
      if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < generalPIDcuts.confNsigmaCombinedProton) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // confNsigmaTPCTOFKaon -> are we doing TPC TOF PID for Kaons? (boolean)
    // confNsigmaTPCKaon -> TPC Kaon Sigma for momentum < 0.5
    // confNsigmaCombinedKaon -> TPC and TOF Kaon Sigma (combined) for momentum > 0.5
    if (mom < 0.3) { // 0.0-0.3
      if (std::abs(nsigmaTPCK) < 3.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (std::abs(nsigmaTPCK) < 2.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (std::abs(nsigmaTPCK) < 1.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((std::abs(nsigmaTOFK) < 3.0) && (std::abs(nsigmaTPCK) < 3.0)) {
        {
          return true;
        }
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((std::abs(nsigmaTOFK) < 2.0) && (std::abs(nsigmaTPCK) < 3.0)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // confNsigmaTPCPion -> TPC Kaon Sigma for momentum < 0.5
    // confNsigmaCombinedPion -> TPC and TOF Pion Sigma (combined) for momentum > 0.5
    if (true) {
      if (mom < 0.5) {
        if (std::abs(nsigmaTPCPi) < generalPIDcuts.confNsigmaTPCPion) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
        if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < generalPIDcuts.confNsigmaCombinedPion) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  bool IsParticleNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    switch (trackonefilter.confPDGCodePartOne) {
      case 2212:  // Proton
      case -2212: // anty Proton
        return IsProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
        break;
      case 211:  // Pion
      case -211: // Pion-
      case 111:  // Pion 0
        return IsPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
        break;
      case 321:  // Kaon+
      case -321: // Kaon-
        return IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
        break;
      default:
        return false;
    }
    return false;
  }

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHisto.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarBins, trackonefilter.confIsMC, trackonefilter.confPDGCodePartOne.value, true);
    vPIDPartOne = confPIDPartOne.value;
    kNsigma = confTrkPIDnSigmaMax.value;
  }

  /// Porduce QA plots for sigle track selection in FemtoUniverse framework
  template <bool isMC, typename PartitionType>
  void FillDebugHistos(o2::aod::FdCollision& col, PartitionType& groupPartsOne)
  {
    eventHisto.fillQA(col);
    for (auto& part : groupPartsOne) {
      // if (part.p() > confCutTable->get("MaxP") || part.pt() > confCutTable->get("MaxPt")) {
      //   continue;
      // }
      // if (!isFullPIDSelected(part.pidCut(), part.p(), confCutTable->get("PIDthr"), vPIDPartOne, confNspecies, kNsigma, confCutTable->get("nSigmaTPC"), confCutTable->get("nSigmaTPCTOF"))) {
      //   continue;
      // }
      if (trackonefilter.confIsTrackIdentified) {
        if (!IsParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
          continue;
        }
      }
      trackHisto.fillQA<isMC, true>(part);
    }
  }

  /// process function when runnning over data/ Monte Carlo reconstructed only
  /// \param col subscribe to FemtoUniverseCollision table
  /// \param parts subscribe to FemtoUniverseParticles table
  void processData(o2::aod::FdCollision& col, FemtoFullParticles&)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    FillDebugHistos<false>(col, groupPartsOne);
  }
  PROCESS_SWITCH(FemtoUniverseDebugTrack, processData, "Enable Debug processing for Monte Carlo", true);

  /// process function when runnning over Monte Carlo with MC truth enabled

  /// \param col subscribe to FemtoUniverseCollision table
  /// \param parts subscribe to the joined table of FemtoUniverseParticles and FemtoUniverseMCLabels table
  /// \param FemtoDramMCParticles subscribe to the table containing the Monte Carlo Truth information
  void processMC(o2::aod::FdCollision& col, FemtoFullParticlesMC&, o2::aod::FdMCParticles&)
  {
    auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    FillDebugHistos<true>(col, groupPartsOne);
  }
  PROCESS_SWITCH(FemtoUniverseDebugTrack, processMC, "Enable Debug processing for Monte Carlo", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniverseDebugTrack>(cfgc),
  };
  return workflow;
}
