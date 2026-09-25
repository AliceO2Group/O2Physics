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
/// \file phiStrangeCorrelator.cxx
/// \brief Table producer for Phi-Strangeness correlation analysis
/// \author Stefano Cannito (stefano.cannito@cern.ch)

// o2-linter: disable=name/workflow-file (File contains multiple necessary structs)

#include "PWGLF/DataModel/LFPhiStrangeCorrelationTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <THn.h>
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PhiMesonCandProducer {
  // Produce the table with the phi candidates information
  Produces<aod::PhimesonCandidatesData> phimesonCandidatesData;
  Produces<aod::PhimesonCandidatesMcReco> phimesonCandidatesMcReco;
  Produces<aod::PhimesonCandidatesMcGen> phimesonCandidatesMcGen;

  HistogramRegistry histos{"phimesonCandidates", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables for phi's daughter tracks selection
  struct : ConfigurableGroup {
    Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0f, "Cut on charge"};
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
    Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
    Configurable<float> cMinKaonPtcut{"cMinKaonPtcut", 0.15f, "Track minimum pt cut"};
    Configurable<float> etaMax{"etaMax", 0.8f, "eta max"};
    Configurable<float> pTToUseTOF{"pTToUseTOF", 0.5f, "pT above which use TOF"};
    Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0f, "Track DCAz cut to PV Maximum"};
    Configurable<std::vector<float>> cMaxDCArToPVPhi{"cMaxDCArToPVPhi", {0.004f, 0.013f, 1.0f}, "Track DCAr cut to PV for Phi"};

    Configurable<bool> cfgIsDCAzParameterized{"cfgIsDCAzParameterized", false, "IsDCAzParameterized"};

    Configurable<float> nSigmaCutTPCKa{"nSigmaCutTPCKa", 3.0f, "Value of the TPC Nsigma cut for Kaons"};
    Configurable<float> nSigmaCutCombinedKa{"nSigmaCutCombinedKa", 3.0f, "Value of the TPC and TOF Nsigma cut for Kaons"};

    Configurable<int> minTPCnClsFound{"minTPCnClsFound", 70, "min number of found TPC clusters"};
  } trackConfigs;

  // Configurables on phi selection
  struct : ConfigurableGroup {
    Configurable<float> maxMPhi{"maxMPhi", 1.2f, "Maximum mass for Phi candidates"};
    Configurable<float> minPhiPt{"minPhiPt", 0.4f, "Minimum pT for Phi candidates"};
    Configurable<float> cfgYAcceptance{"cfgYAcceptance", 0.5f, "Rapidity acceptance"};
  } phiConfigs;

  Configurable<bool> trueGenPhi{"trueGenPhi", false, "Phi trigger in MC Collisions: false - gen K+K- pair, true - gen phi meson"};

  // Configurables on phi pT bins
  Configurable<std::vector<double>> binspTPhi{"binspTPhi", {0.4, 0.8, 1.4, 2.0, 2.8, 4.0, 6.0, 10.0}, "pT bin limits for Phi"};

  // Filter on default selected collisions
  // Filter collisionFilter = aod::lf_selection_default_collision::defaultSel == true;
  // Filter collisionFilter = aod::lf_selection_event::defaultSel == true;

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::PVMults>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;

  /*using SelCollisions = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::PVMults, aod::PhiStrangeDefEvtSelDataLike>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels, aod::PhiStrangeDefEvtSelMcGen>;

  using FilteredSelCollisions = soa::Filtered<SelCollisions>;
  using FilteredSimCollisions = soa::Filtered<SimCollisions>;*/

  // Defining the type of the phi's daughter tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;

  Partition<FullTracks> posTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FullTracks> negTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  Partition<FullMCTracks> posMCTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FullMCTracks> negMCTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  // Manual slicing
  Preslice<aod::Tracks> trackPerCollision = aod::track::collisionId;
  SliceCache cache;

  // Constants
  double massKa = o2::constants::physics::MassKPlus;

  void init(InitContext&)
  {
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};
    AxisSpec binnedpTPhiAxis{(std::vector<double>)binspTPhi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massPhiAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};

    histos.add("h3PhiCandidatesMass", "Phi candidate invariant mass", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis});

    // Phi's daughter tracks information
    /*histos.add("hEta", "Eta of Kaon candidates", kTH1F, {{100, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});
    histos.add("h2DauTracksPhiDCAxy", "DCAxy distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    histos.add("h2DauTracksPhiDCAz", "DCAz distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});*/
  }

  // Topological track selection
  template <typename T>
  bool selectionTrackResonance(const T& track)
  {
    if (trackConfigs.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (trackConfigs.cfgPVContributor && !track.isPVContributor()) {
      return false;
    }

    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound) {
      return false;
    }

    if (track.pt() < trackConfigs.cMinKaonPtcut) {
      return false;
    }

    if (std::abs(track.dcaXY()) > trackConfigs.cMaxDCArToPVPhi->at(0) + (trackConfigs.cMaxDCArToPVPhi->at(1) / std::pow(track.pt(), trackConfigs.cMaxDCArToPVPhi->at(2)))) {
      return false;
    }
    if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVcut) {
      return false;
    }

    return true;
  }

  // PIDKaon track selection
  template <typename T>
  bool selectionPIDKaonpTdependent(const T& track)
  {
    if (track.pt() < trackConfigs.pTToUseTOF && std::abs(track.tpcNSigmaKa()) < trackConfigs.nSigmaCutTPCKa) {
      return true;
    }
    if (track.pt() >= trackConfigs.pTToUseTOF && track.hasTOF() && (std::pow(track.tofNSigmaKa(), 2) + std::pow(track.tpcNSigmaKa(), 2)) < std::pow(trackConfigs.nSigmaCutCombinedKa, 2)) {
      return true;
    }

    return false;
  }

  // Reconstruct the Phi candidate
  template <typename T>
  ROOT::Math::PxPyPzMVector recMother(const T& track1, const T& track2, float massdaughter1, float massdaughter2)
  {
    ROOT::Math::PxPyPzMVector daughter1(track1.px(), track1.py(), track1.pz(), massdaughter1); // set the daughter1 4-momentum
    ROOT::Math::PxPyPzMVector daughter2(track2.px(), track2.py(), track2.pz(), massdaughter2); // set the daughter2 4-momentum
    ROOT::Math::PxPyPzMVector mother = daughter1 + daughter2;                                  // calculate the mother 4-momentum

    return mother;
  }

  template <bool isMC, typename T>
  void processPhiCandidates(const T& collision, std::optional<std::reference_wrapper<const aod::McParticles>> mcParticlesOpt = std::nullopt)
  {
    // Compile-time selection of the track partition to use based on the type of the analysis (Data or MCReco)
    auto posThisColl = [&]() {
      if constexpr (isMC) {
        return posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      } else {
        return posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      }
    }();

    auto negThisColl = [&]() {
      if constexpr (isMC) {
        return negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      } else {
        return negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      }
    }();

    // Loop over positive tracks
    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1)) {
        continue;
      }

      if constexpr (isMC) {
        if (!track1.has_mcParticle()) {
          continue;
        }
        if (!mcParticlesOpt) {
          continue;
        }

        const auto& mcParticles = mcParticlesOpt->get();
        const auto track1McParticle = mcParticles.rawIteratorAt(track1.mcParticleId());
        if (track1McParticle.pdgCode() != PDG_t::kKPlus || !track1McParticle.isPhysicalPrimary()) {
          continue;
        }
      }

      // Loop over negative tracks
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2)) {
          continue;
        }

        if constexpr (isMC) {
          if (!track2.has_mcParticle()) {
            continue;
          }
          if (!mcParticlesOpt) {
            continue;
          }

          const auto& mcParticles = mcParticlesOpt->get();
          const auto track2McParticle = mcParticles.rawIteratorAt(track2.mcParticleId());
          if (track2McParticle.pdgCode() != PDG_t::kKMinus || !track2McParticle.isPhysicalPrimary()) {
            continue;
          }
        }

        // Kinematic reconstruction (common for data and MC)
        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);

        if (recPhi.Pt() < phiConfigs.minPhiPt) {
          continue;
        }
        if (recPhi.M() > phiConfigs.maxMPhi) {
          continue;
        }
        if (std::abs(recPhi.Rapidity()) > phiConfigs.cfgYAcceptance) {
          continue;
        }

        // PDG check and MC truth association for MCReco analysis
        if constexpr (isMC) {
          if (!mcParticlesOpt) {
            continue;
          }

          const auto& mcParticles = mcParticlesOpt->get();

          const auto track1McParticle = mcParticles.rawIteratorAt(track1.mcParticleId());
          const auto track2McParticle = mcParticles.rawIteratorAt(track2.mcParticleId());

          const auto track1mcPartMotherIndexes = track1McParticle.mothersIds();
          const auto track2mcPartMotherIndexes = track2McParticle.mothersIds();

          auto genPhiMaybe = [&]() -> std::optional<aod::McParticles::iterator> {
            for (const auto& mother1Index : track1mcPartMotherIndexes) {
              for (const auto& mother2Index : track2mcPartMotherIndexes) {
                if (mother1Index != mother2Index) {
                  continue;
                }
                auto motherMcParticle = mcParticles.rawIteratorAt(mother1Index);
                if (std::abs(motherMcParticle.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
                  return motherMcParticle;
                }
              }
            }
            return std::nullopt;
          }();

          if (!genPhiMaybe) {
            continue;
          }
          const auto& genPhi = *genPhiMaybe;

          phimesonCandidatesMcReco(collision.globalIndex(), recPhi.M(), genPhi.pt(), genPhi.y(), genPhi.phi());
        } else {
          histos.fill(HIST("h3PhiCandidatesMass"), collision.centFT0M(), recPhi.Pt(), recPhi.M());
          phimesonCandidatesData(collision.globalIndex(), recPhi.M(), recPhi.Pt(), recPhi.Rapidity(), recPhi.Phi());
        }
      }
    }
  }

  void processData(SelCollisions::iterator const& collision, FullTracks const&)
  {
    processPhiCandidates<false>(collision);
  }

  PROCESS_SWITCH(PhiMesonCandProducer, processData, "Process function to select Phi meson candidates in Data or in McReco (w/o McTruth) analysis", true);

  void processMCReco(SimCollisions::iterator const& collision, FullMCTracks const&, aod::McParticles const& mcParticles)
  {
    processPhiCandidates<true>(collision, mcParticles);
  }

  PROCESS_SWITCH(PhiMesonCandProducer, processMCReco, "Process function to select Phi meson candidates in MCReco w MC truth", false);

  /*void processFilteredData(FilteredSelCollisions::iterator const& collision, FullTracks const&)
  {
    processPhiCandidates<false>(collision);
  }

  PROCESS_SWITCH(PhiMesonCandProducer, processFilteredData, "Process function to select Phi meson candidates in filtered Data", true);

  void processFilteredMCReco(FilteredSimCollisions::iterator const& collision, FullMCTracks const&, aod::McParticles const& mcParticles)
  {
    processPhiCandidates<true>(collision, mcParticles);
  }

  PROCESS_SWITCH(PhiMesonCandProducer, processFilteredMCReco, "Process function to select Phi meson candidates in filtered MCReco w MC truth", false);*/

  void processMCGen(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    if (trueGenPhi) {
      for (const auto& mcParticle : mcParticles) {
        if (std::abs(mcParticle.pdgCode()) != o2::constants::physics::Pdg::kPhi) {
          continue;
        }
        if (mcParticle.pt() < phiConfigs.minPhiPt) {
          continue;
        }
        if (std::abs(mcParticle.y()) > phiConfigs.cfgYAcceptance) {
          continue;
        }

        phimesonCandidatesMcGen(mcCollision.globalIndex(), 0, mcParticle.pt(), mcParticle.y(), mcParticle.phi());
      }
    } else {
      for (const auto& mcParticle1 : mcParticles) {
        if (!mcParticle1.isPhysicalPrimary() || std::abs(mcParticle1.eta()) > trackConfigs.etaMax) {
          continue;
        }

        for (const auto& mcParticle2 : mcParticles) {
          if (!mcParticle2.isPhysicalPrimary() || std::abs(mcParticle2.eta()) > trackConfigs.etaMax) {
            continue;
          }

          if ((mcParticle1.pdgCode() != PDG_t::kKPlus || mcParticle2.pdgCode() != PDG_t::kKMinus) &&
              (mcParticle1.pdgCode() != PDG_t::kKMinus || mcParticle2.pdgCode() != PDG_t::kKPlus)) {
            continue;
          }

          ROOT::Math::PxPyPzMVector genKPair = recMother(mcParticle1, mcParticle2, massKa, massKa);

          if (genKPair.Pt() < phiConfigs.minPhiPt) {
            continue;
          }
          if (genKPair.M() > phiConfigs.maxMPhi) {
            continue;
          }
          if (std::abs(genKPair.Rapidity()) > phiConfigs.cfgYAcceptance) {
            continue;
          }

          phimesonCandidatesMcGen(mcCollision.globalIndex(), genKPair.M(), genKPair.Pt(), genKPair.Rapidity(), genKPair.Phi());
        }
      }
    }
  }

  PROCESS_SWITCH(PhiMesonCandProducer, processMCGen, "Process function to select Phi meson candidates in MCGen", false);
};

struct V0ReducedCandProducer {
  enum V0Type { kK0S = 0,
                kLambda,
                kAntiLambda };

  // Produce the table with the V0 candidates information
  Produces<aod::K0sReducedCandidatesData> k0sReducedCandidatesData;
  Produces<aod::K0sReducedCandidatesMcReco> k0sReducedCandidatesMcReco;

  Produces<aod::LambdaReducedCandidatesData> lambdaReducedCandidatesData;
  Produces<aod::LambdaReducedCandidatesMcReco> lambdaReducedCandidatesMcReco;

  Produces<aod::AntiLambdaReducedCandidatesData> antilambdaReducedCandidatesData;
  Produces<aod::AntiLambdaReducedCandidatesMcReco> antilambdaReducedCandidatesMcReco;

  HistogramRegistry histos{"k0sReducedCandidates", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables for tracks selection
  struct : ConfigurableGroup {
    Configurable<float> etaMax{"etaMax", 0.8f, "eta max"};
    Configurable<float> nSigmaCutTPCSecPion{"nSigmaCutTPCSecPion", 4.0f, "Value of the TPC Nsigma cut for secondary Pions"};
    Configurable<float> nSigmaCutTPCSecProton{"nSigmaCutTPCSecProton", 4.0f, "Value of the TPC Nsigma cut for secondary Protons"};

    Configurable<int> minTPCnClsFound{"minTPCnClsFound", 70, "min number of found TPC clusters"};
    Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70, "min number of TPC crossed rows"};
    Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  } trackConfigs;

  // Configurables for V0 selection
  struct : ConfigurableGroup {
    Configurable<float> decayRadius{"decayRadius", 0.5f, "V0 decay radius"};
    Configurable<float> dcaV0Dau{"dcaV0Dau", 1.0f, "DCA V0 Daughters"};
    Configurable<float> dcaPosToPV{"dcaPosToPV", 0.1f, "DCA Pos To PV"};
    Configurable<float> dcaNegToPV{"dcaNegToPV", 0.1f, "DCA Neg To PV"};

    Configurable<bool> furtherV0Selection{"furtherV0Selection", false, "Further V0 selection"};

    Configurable<float> yAcceptance{"yAcceptance", 0.5f, "Rapidity acceptance"};

    struct : ConfigurableGroup {
      Configurable<float> cosPAK0S{"cosPAK0S", 0.98f, "K0S CosPA"};
      Configurable<float> ctauK0S{"ctauK0S", 20.0f, "C tau K0s(cm)"};
      Configurable<float> rejMassK0S{"rejMassK0S", 0.005f, "V0 rej K0s"};

      Configurable<float> paramArmenterosCutK0S{"paramArmenterosCutK0S", 0.2f, "Parameter for Armenteros Cut for K0S"};

      Configurable<float> minPtK0S{"minPtK0S", 0.1f, "K0S min pt"};
      Configurable<std::vector<double>> binspTK0S{"binspTK0S", {0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for K0S"};
      Configurable<std::pair<float, float>> rangeMK0sSignal{"rangeMK0sSignal", {0.47f, 0.53f}, "K0S mass range for signal extraction"};
    } k0sConfigs;

    struct : ConfigurableGroup {
      Configurable<float> cosPALambda{"cosPALambda", 0.98f, "Lambda CosPA"};
      Configurable<float> ctauLambda{"ctauLambda", 7.89f, "C tau Lambda(cm)"};
      Configurable<float> rejMassLambda{"rejMassLambda", 0.005f, "V0 rej Lambda"};

      Configurable<float> minPtLambda{"minPtLambda", 0.1f, "Lambda min pt"};
      Configurable<std::vector<double>> binspTLambda{"binspTLambda", {0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Lambda"};
      Configurable<std::pair<float, float>> rangeMLambdaSignal{"rangeMLambdaSignal", {1.11f, 1.12f}, "Lambda mass range for signal extraction"};
    } lambdaConfigs;
  } v0Configs;

  Configurable<std::vector<double>> binspTXi{"binspTXi", {0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Xi"};

  // Constants
  static constexpr double MassK0S = o2::constants::physics::MassK0Short;
  static constexpr double MassLambda = o2::constants::physics::MassLambda0;

  // Filter on default selected collisions
  // Filter collisionFilter = aod::lf_selection_default_collision::defaultSel == true;
  Filter collisionFilter = aod::lf_selection_event::defaultSel == true;

  // Defining filters on V0s (cannot filter on dynamic columns)
  Filter v0PreFilter = (nabs(aod::v0data::dcapostopv) > v0Configs.dcaPosToPV &&
                        nabs(aod::v0data::dcanegtopv) > v0Configs.dcaNegToPV &&
                        aod::v0data::dcaV0daughters < v0Configs.dcaV0Dau);

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::PVMults, aod::PhiStrangeEvtSelDataLike>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;

  using FilteredSelCollisions = soa::Filtered<SelCollisions>;
  using FilteredSimCollisions = soa::Filtered<SimCollisions>;

  // Defining the type of the V0s and corresponding daughter tracks for data and MC
  using FullV0s = soa::Filtered<aod::V0Datas>;
  using FullMCV0s = soa::Filtered<soa::Join<aod::V0Datas, aod::McV0Labels>>;

  using V0DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using V0DauMCTracks = soa::Join<V0DauTracks, aod::McTrackLabels>;

  void init(InitContext&)
  {
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};
    AxisSpec binnedpTK0SAxis{(std::vector<double>)v0Configs.k0sConfigs.binspTK0S, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTLambdaAxis{(std::vector<double>)v0Configs.lambdaConfigs.binspTLambda, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTXiAxis{(std::vector<double>)binspTXi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massK0SAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec massLambdaAxis = {200, 1.08f, 1.15f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};

    histos.add("h3K0sCandidatesMass", "K^{0}_{S} candidate invariant mass", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massK0SAxis});
    histos.add("h3LambdaCandidatesMass", "#Lambda candidate invariant mass", kTH3F, {binnedmultAxis, binnedpTLambdaAxis, massLambdaAxis});
    histos.add("h3AntiLambdaCandidatesMass", "#bar{#Lambda} candidate invariant mass", kTH3F, {binnedmultAxis, binnedpTLambdaAxis, massLambdaAxis});

    histos.add("h2LambdaFromXi", "#Lambda from #Xi", kTH2F, {binnedpTLambdaAxis, binnedpTXiAxis});
    histos.add("h2AntiLambdaFromXi", "#bar{#Lambda} from #Xi", kTH2F, {binnedpTLambdaAxis, binnedpTXiAxis});
  }

  template <V0Type v0Type>
  constexpr int getPdgCode()
  {
    if constexpr (v0Type == kK0S) {
      return PDG_t::kK0Short;
    } else if constexpr (v0Type == kLambda) {
      return PDG_t::kLambda0;
    } else if constexpr (v0Type == kAntiLambda) {
      return PDG_t::kLambda0Bar;
    } else {
      static_assert(v0Type == kK0S || v0Type == kLambda || v0Type == kAntiLambda, "Unsupported particle type in getPdgCode");
      return 0;
    }
  }

  // Single track selection for strangeness sector
  template <typename T>
  bool selectionTrackStrangeness(const T& track)
  {
    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < trackConfigs.minNCrossedRowsTPC) {
      return false;
    }
    if (track.tpcChi2NCl() > trackConfigs.maxChi2TPC) {
      return false;
    }

    if (std::abs(track.eta()) > trackConfigs.etaMax) {
      return false;
    }
    return true;
  }

  // V0 selection
  template <bool isMC, V0Type v0Type, typename T1, typename T2>
  bool selectionV0(const T1& v0, const T2& collision)
  {
    using V0DauTrackType = std::conditional_t<isMC, V0DauMCTracks, V0DauTracks>;

    const auto& posDaughterTrack = v0.template posTrack_as<V0DauTrackType>();
    const auto& negDaughterTrack = v0.template negTrack_as<V0DauTrackType>();

    if (!selectionTrackStrangeness(posDaughterTrack) || !selectionTrackStrangeness(negDaughterTrack)) {
      return false;
    }

    if constexpr (!isMC) {
      if constexpr (v0Type == kK0S) {
        if (std::abs(posDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCSecPion) {
          return false;
        }
        if (std::abs(negDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCSecPion) {
          return false;
        }
      } else if constexpr (v0Type == kLambda) {
        if (std::abs(posDaughterTrack.tpcNSigmaPr()) > trackConfigs.nSigmaCutTPCSecProton) {
          return false;
        }
        if (std::abs(negDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCSecPion) {
          return false;
        }
      } else if constexpr (v0Type == kAntiLambda) {
        if (std::abs(posDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCSecPion) {
          return false;
        }
        if (std::abs(negDaughterTrack.tpcNSigmaPr()) > trackConfigs.nSigmaCutTPCSecProton) {
          return false;
        }
      }
    }

    if (v0.v0radius() < v0Configs.decayRadius) {
      return false;
    }
    if (std::abs(v0.rapidity(v0Type)) > v0Configs.yAcceptance) {
      return false;
    }

    const float cutCosPA = (v0Type == kK0S) ? v0Configs.k0sConfigs.cosPAK0S : v0Configs.lambdaConfigs.cosPALambda;
    const float cutMinPt = (v0Type == kK0S) ? v0Configs.k0sConfigs.minPtK0S : v0Configs.lambdaConfigs.minPtLambda;
    const float cutCtau = (v0Type == kK0S) ? v0Configs.k0sConfigs.ctauK0S : v0Configs.lambdaConfigs.ctauLambda;

    if (v0.v0cosPA() < cutCosPA) {
      return false;
    }
    if (v0.pt() < cutMinPt) {
      return false;
    }

    if (v0Configs.furtherV0Selection) {
      constexpr float assumedMass = v0Type == kK0S ? MassK0S : MassLambda;
      if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * assumedMass > cutCtau) {
        return false;
      }

      if constexpr (v0Type == kK0S) {
        if (std::abs(v0.mLambda() - MassLambda) < v0Configs.k0sConfigs.rejMassK0S ||
            std::abs(v0.mAntiLambda() - MassLambda) < v0Configs.k0sConfigs.rejMassK0S) {
          return false;
        }
        if (v0.qtarm() < (v0Configs.k0sConfigs.paramArmenterosCutK0S * std::abs(v0.alpha()))) {
          return false;
        }
      } else {
        if (std::abs(v0.mK0Short() - MassK0S) < v0Configs.lambdaConfigs.rejMassLambda) {
          return false;
        }
      }
    }

    return true;
  }

  template <bool isMC, V0Type v0Type, typename TV0, typename TColl>
  void processV0Type(TColl const& collision, TV0 const& v0, std::optional<std::reference_wrapper<const aod::McParticles>> mcParticlesOpt = std::nullopt)
  {
    // Cut on V0 dynamic columns
    if (!selectionV0<isMC, v0Type>(v0, collision)) {
      return;
    }

    if constexpr (isMC) {
      if (!v0.has_mcParticle() || !mcParticlesOpt) {
        return;
      }

      const auto& mcParticles = mcParticlesOpt->get();
      const auto& v0McParticle = mcParticles.rawIteratorAt(v0.mcParticleId());
      // if (v0McParticle.pdgCode() != getPdgCode<v0Type>() || !v0McParticle.isPhysicalPrimary()) {
      // return;
      //}
      if (v0McParticle.pdgCode() != getPdgCode<v0Type>()) {
        return;
      }

      if constexpr ((v0Type == kLambda || v0Type == kAntiLambda)) {
        if (v0.has_mcMotherParticle()) {
          auto v0McMotherParticle = mcParticles.rawIteratorAt(v0.mcMotherParticleId());
          if (std::abs(v0McMotherParticle.y()) <= v0Configs.yAcceptance) {
            if constexpr (v0Type == kLambda) {
              if (v0McMotherParticle.pdgCode() == PDG_t::kXiMinus || v0McMotherParticle.pdgCode() == o2::constants::physics::Pdg::kXi0) {
                histos.fill(HIST("h2LambdaFromXi"), v0McParticle.pt(), v0McMotherParticle.pt());
              }
            } else if constexpr (v0Type == kAntiLambda) {
              if (v0McMotherParticle.pdgCode() == PDG_t::kXiPlusBar || v0McMotherParticle.pdgCode() == o2::constants::physics::Pdg::kXi0) {
                histos.fill(HIST("h2AntiLambdaFromXi"), v0McParticle.pt(), v0McMotherParticle.pt());
              }
            }
          }
        }
      }

      if (!v0McParticle.isPhysicalPrimary()) {
        return;
      }

      if constexpr (v0Type == kK0S) {
        k0sReducedCandidatesMcReco(collision.globalIndex(), v0.mK0Short(), v0McParticle.pt(), v0McParticle.y(), v0McParticle.phi());
      } else if constexpr (v0Type == kLambda) {
        lambdaReducedCandidatesMcReco(collision.globalIndex(), v0.mLambda(), v0McParticle.pt(), v0McParticle.y(), v0McParticle.phi());
      } else if constexpr (v0Type == kAntiLambda) {
        antilambdaReducedCandidatesMcReco(collision.globalIndex(), v0.mAntiLambda(), v0McParticle.pt(), v0McParticle.y(), v0McParticle.phi());
      }
    } else {
      if constexpr (v0Type == kK0S) {
        histos.fill(HIST("h3K0sCandidatesMass"), collision.centFT0M(), v0.pt(), v0.mK0Short());
        k0sReducedCandidatesData(collision.globalIndex(), v0.mK0Short(), v0.pt(), v0.yK0Short(), v0.phi());
      } else if constexpr (v0Type == kLambda) {
        histos.fill(HIST("h3LambdaCandidatesMass"), collision.centFT0M(), v0.pt(), v0.mLambda());
        lambdaReducedCandidatesData(collision.globalIndex(), v0.mLambda(), v0.pt(), v0.yLambda(), v0.phi());
      } else if constexpr (v0Type == kAntiLambda) {
        histos.fill(HIST("h3AntiLambdaCandidatesMass"), collision.centFT0M(), v0.pt(), v0.mAntiLambda());
        antilambdaReducedCandidatesData(collision.globalIndex(), v0.mAntiLambda(), v0.pt(), v0.yLambda(), v0.phi());
        // histos.fill(HIST("h3LambdaCandidatesMass"), collision.centFT0M(), v0.pt(), v0.mAntiLambda());
        // lambdaReducedCandidatesData(collision.globalIndex(), v0.mAntiLambda(), v0.pt(), v0.yLambda(), v0.phi());
      }
    }
  }

  template <bool isMC, V0Type... Types, typename TColl, typename TV0>
  void evaluateAllHypotheses(TColl const& collision, TV0 const& v0, std::optional<std::reference_wrapper<const aod::McParticles>> mcParticlesOpt = std::nullopt)
  {
    (processV0Type<isMC, Types>(collision, v0, mcParticlesOpt), ...);
  }

  // void processData(SelCollisions::iterator const& collision, FullV0s const& V0s, V0DauTracks const&)
  void processData(FilteredSelCollisions::iterator const& collision, FullV0s const& V0s, V0DauTracks const&)
  {
    for (const auto& v0 : V0s) {
      evaluateAllHypotheses<false, kK0S, kLambda, kAntiLambda>(collision, v0);
    }
  }

  PROCESS_SWITCH(V0ReducedCandProducer, processData, "Process function to select reduced K0S candidates in Data or in McReco (w/o McTruth) analysis", true);

  // void processMCReco(SimCollisions::iterator const& collision, FullMCV0s const& V0s, V0DauMCTracks const&, aod::McParticles const& mcParticles)
  void processMCReco(FilteredSimCollisions::iterator const& collision, FullMCV0s const& V0s, V0DauMCTracks const&, aod::McParticles const& mcParticles)
  {
    for (const auto& v0 : V0s) {
      evaluateAllHypotheses<true, kK0S, kLambda, kAntiLambda>(collision, v0, mcParticles);
    }
  }

  PROCESS_SWITCH(V0ReducedCandProducer, processMCReco, "Process function to select reduced K0s candidates in MCReco w MC truth", false);
};

struct CascadeReducedCandProducer {
  enum CascadeType { kXi = 0,
                     kAntiXi,
                     kOmega,
                     kAntiOmega };

  // Produce the table with the Xi candidates information
  Produces<aod::XiReducedCandidatesData> xiReducedCandidatesData;
  Produces<aod::XiReducedCandidatesMcReco> xiReducedCandidatesMcReco;

  Produces<aod::OmegaReducedCandidatesData> omegaReducedCandidatesData;
  Produces<aod::OmegaReducedCandidatesMcReco> omegaReducedCandidatesMcReco;

  HistogramRegistry histos{"xiReducedCandidates", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables for tracks selection
  struct : ConfigurableGroup {
    Configurable<float> etaMax{"etaMax", 0.8f, "eta max"};
    Configurable<float> nSigmaCutTPCSecPion{"nSigmaCutTPCSecPion", 4.0f, "Value of the TPC Nsigma cut for secondary Pions"};
    Configurable<float> nSigmaCutTPCSecKaon{"nSigmaCutTPCSecKaon", 4.0f, "Value of the TPC Nsigma cut for secondary Kaons"};
    Configurable<float> nSigmaCutTPCSecProton{"nSigmaCutTPCSecProton", 4.0f, "Value of the TPC Nsigma cut for secondary Protons"};

    Configurable<int> minTPCnClsFound{"minTPCnClsFound", 70, "min number of found TPC clusters"};
    Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70, "min number of TPC crossed rows"};
    Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  } trackConfigs;

  // Configurables for Cascade selection
  struct : ConfigurableGroup {
    Configurable<float> dcaV0Dau{"dcaV0Dau", 1.0f, "DCA V0 Daughters"};
    Configurable<float> dcaCascDau{"dcaCascDau", 1.0f, "DCA Cascade Daughters"};
    Configurable<float> dcaPosToPV{"dcaPosToPV", 0.1f, "DCA Pos To PV"};
    Configurable<float> dcaNegToPV{"dcaNegToPV", 0.1f, "DCA Neg To PV"};
    Configurable<float> dcaBachToPV{"dcaBachToPV", 0.1f, "DCA Bach To PV"};

    Configurable<float> yAcceptance{"yAcceptance", 0.5f, "Rapidity acceptance"};

    struct : ConfigurableGroup {
      Configurable<float> v0CosPAXi{"v0CosPAXi", 0.98f, "Xi V0 CosPA"};
      Configurable<float> cosPAXi{"cosPAXi", 0.98f, "Xi CosPA"};
      Configurable<float> v0RadiusXi{"v0RadiusXi", 0.5f, "V0 decay radius"};
      Configurable<float> radiusXi{"radiusXi", 0.5f, "Xi decay radius"};

      Configurable<float> minPtXi{"minPtXi", 0.8f, "Xi min pt"};
      Configurable<std::vector<double>> binspTXi{"binspTXi", {0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Xi"};
    } xiConfigs;

    struct : ConfigurableGroup {
      Configurable<float> v0CosPAOmega{"v0CosPAOmega", 0.98f, "Omega V0 CosPA"};
      Configurable<float> cosPAOmega{"cosPAOmega", 0.98f, "Omega CosPA"};
      Configurable<float> v0RadiusOmega{"v0RadiusOmega", 0.5f, "V0 decay radius"};
      Configurable<float> radiusOmega{"radiusOmega", 0.5f, "Omega decay radius"};

      Configurable<float> minPtOmega{"minPtOmega", 0.8f, "Omega min pt"};
      Configurable<std::vector<double>> binspTOmega{"binspTOmega", {0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Omega"};
    } omegaConfigs;
  } cascadeConfigs;

  // Constants
  static constexpr double MassXi = o2::constants::physics::MassXiMinus;
  static constexpr double MassOmega = o2::constants::physics::MassOmegaMinus;

  // Filter on default selected collisions
  Filter collisionFilter = aod::lf_selection_event::defaultSel == true;

  // Defining filters on Cascades (cannot filter on dynamic columns)
  Filter cascadePreFilter = (nabs(aod::cascdata::dcapostopv) > cascadeConfigs.dcaPosToPV &&
                             nabs(aod::cascdata::dcanegtopv) > cascadeConfigs.dcaNegToPV &&
                             nabs(aod::cascdata::dcabachtopv) > cascadeConfigs.dcaBachToPV &&
                             aod::cascdata::dcaV0daughters < cascadeConfigs.dcaV0Dau &&
                             aod::cascdata::dcacascdaughters < cascadeConfigs.dcaCascDau);

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::PVMults, aod::PhiStrangeEvtSelDataLike>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;

  using FilteredSelCollisions = soa::Filtered<SelCollisions>;
  using FilteredSimCollisions = soa::Filtered<SimCollisions>;

  // Defining the type of the Cascades and corresponding daughter and bachelor tracks for data and MC
  using FullCascades = soa::Filtered<aod::CascDatas>;
  using FullMCCascades = soa::Filtered<soa::Join<aod::CascDatas, aod::McCascLabels>>;

  using DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
  using DauMCTracks = soa::Join<DauTracks, aod::McTrackLabels>;

  void init(InitContext&)
  {
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};
    AxisSpec binnedpTXiAxis{(std::vector<double>)cascadeConfigs.xiConfigs.binspTXi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTOmegaAxis{(std::vector<double>)cascadeConfigs.omegaConfigs.binspTOmega, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massXiAxis = {200, 1.2f, 1.4f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec massOmegaAxis = {200, 1.6f, 1.8f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};

    histos.add("h3XiCandidatesMass", "#Xi candidate invariant mass", kTH3F, {binnedmultAxis, binnedpTXiAxis, massXiAxis});
    histos.add("h3OmegaCandidatesMass", "#Omega candidate invariant mass", kTH3F, {binnedmultAxis, binnedpTOmegaAxis, massOmegaAxis});
  }

  template <CascadeType cascadeType>
  constexpr int getPdgCode()
  {
    if constexpr (cascadeType == kXi) {
      return PDG_t::kXiMinus;
    } else if constexpr (cascadeType == kAntiXi) {
      return PDG_t::kXiPlusBar;
    } else if constexpr (cascadeType == kOmega) {
      return PDG_t::kOmegaMinus;
    } else if constexpr (cascadeType == kAntiOmega) {
      return PDG_t::kOmegaPlusBar;
    } else {
      static_assert(cascadeType == kXi || cascadeType == kAntiXi || cascadeType == kOmega || cascadeType == kAntiOmega, "Unsupported particle type in getPdgCode");
      return 0;
    }
  }

  // Single track selection for strangeness sector
  template <typename T>
  bool selectionTrackStrangeness(const T& track)
  {
    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < trackConfigs.minNCrossedRowsTPC) {
      return false;
    }
    if (track.tpcChi2NCl() > trackConfigs.maxChi2TPC) {
      return false;
    }

    if (std::abs(track.eta()) > trackConfigs.etaMax) {
      return false;
    }
    return true;
  }

  // Cascade selection
  template <bool isMC, CascadeType cascadeType, typename T1, typename T2>
  bool selectionCascade(const T1& cascade, const T2& collision)
  {
    if constexpr (cascadeType == kXi || cascadeType == kOmega) {
      if (cascade.sign() > 0) {
        return false; // Matter must have negative charge
      }
    } else {
      if (cascade.sign() < 0) {
        return false; // Anti-matter must have positive charge
      }
    }

    using DauTrackType = std::conditional_t<isMC, DauMCTracks, DauTracks>;

    const auto& posDaughterTrack = cascade.template posTrack_as<DauTrackType>();
    const auto& negDaughterTrack = cascade.template negTrack_as<DauTrackType>();
    const auto& bachDaughterTrack = cascade.template bachelor_as<DauTrackType>();

    if (!selectionTrackStrangeness(posDaughterTrack) || !selectionTrackStrangeness(negDaughterTrack) || !selectionTrackStrangeness(bachDaughterTrack)) {
      return false;
    }

    if constexpr (!isMC) {
      if constexpr (cascadeType == kXi || cascadeType == kOmega) {
        if (std::abs(posDaughterTrack.tpcNSigmaPr()) > trackConfigs.nSigmaCutTPCSecProton) {
          return false;
        }
        if (std::abs(negDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCSecPion) {
          return false;
        }
      } else {
        if (std::abs(posDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCSecPion) {
          return false;
        }
        if (std::abs(negDaughterTrack.tpcNSigmaPr()) > trackConfigs.nSigmaCutTPCSecProton) {
          return false;
        }
      }

      if constexpr (cascadeType == kXi || cascadeType == kAntiXi) {
        if (std::abs(bachDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCSecPion) {
          return false;
        }
      } else if constexpr (cascadeType == kOmega || cascadeType == kAntiOmega) {
        if (std::abs(bachDaughterTrack.tpcNSigmaKa()) > trackConfigs.nSigmaCutTPCSecKaon) {
          return false;
        }
      }
    }

    if (std::abs(cascade.rapidity(cascadeType)) > cascadeConfigs.yAcceptance) {
      return false;
    }

    const float cutV0CosPA = (cascadeType == kXi || cascadeType == kAntiXi) ? cascadeConfigs.xiConfigs.v0CosPAXi : cascadeConfigs.omegaConfigs.v0CosPAOmega;
    const float cutCosPA = (cascadeType == kXi || cascadeType == kAntiXi) ? cascadeConfigs.xiConfigs.cosPAXi : cascadeConfigs.omegaConfigs.cosPAOmega;
    const float cutV0Radius = (cascadeType == kXi || cascadeType == kAntiXi) ? cascadeConfigs.xiConfigs.v0RadiusXi : cascadeConfigs.omegaConfigs.v0RadiusOmega;
    const float cutRadius = (cascadeType == kXi || cascadeType == kAntiXi) ? cascadeConfigs.xiConfigs.radiusXi : cascadeConfigs.omegaConfigs.radiusOmega;
    const float cutMinPt = (cascadeType == kXi || cascadeType == kAntiXi) ? cascadeConfigs.xiConfigs.minPtXi : cascadeConfigs.omegaConfigs.minPtOmega;

    const auto& pvx = collision.posX();
    const auto& pvy = collision.posY();
    const auto& pvz = collision.posZ();

    if (cascade.v0cosPA(pvx, pvy, pvz) < cutV0CosPA) {
      return false;
    }
    if (cascade.casccosPA(pvx, pvy, pvz) < cutCosPA) {
      return false;
    }
    if (cascade.v0radius() < cutV0Radius) {
      return false;
    }
    if (cascade.cascradius() < cutRadius) {
      return false;
    }
    if (cascade.pt() < cutMinPt) {
      return false;
    }

    return true;
  }

  template <bool isMC, CascadeType cascadeType, typename TCascade, typename TColl>
  void processCascadeType(TColl const& collision, TCascade const& cascade, std::optional<std::reference_wrapper<const aod::McParticles>> mcParticlesOpt = std::nullopt)
  {
    // Cut on Cascade dynamic columns
    if (!selectionCascade<isMC, cascadeType>(cascade, collision)) {
      return;
    }

    if constexpr (isMC) {
      if (!cascade.has_mcParticle() || !mcParticlesOpt) {
        return;
      }

      const auto& mcParticles = mcParticlesOpt->get();
      const auto& cascadeMcParticle = mcParticles.rawIteratorAt(cascade.mcParticleId());
      if (cascadeMcParticle.pdgCode() != getPdgCode<cascadeType>() || !cascadeMcParticle.isPhysicalPrimary()) {
        return;
      }

      if constexpr (cascadeType == kXi || cascadeType == kAntiXi) {
        xiReducedCandidatesMcReco(collision.globalIndex(), cascade.mXi(), cascadeMcParticle.pt(), cascadeMcParticle.y(), cascadeMcParticle.phi());
      } else if constexpr (cascadeType == kOmega || cascadeType == kAntiOmega) {
        omegaReducedCandidatesMcReco(collision.globalIndex(), cascade.mOmega(), cascadeMcParticle.pt(), cascadeMcParticle.y(), cascadeMcParticle.phi());
      }
    } else {
      if constexpr (cascadeType == kXi || cascadeType == kAntiXi) {
        histos.fill(HIST("h3XiCandidatesMass"), collision.centFT0M(), cascade.pt(), cascade.mXi());
        xiReducedCandidatesData(collision.globalIndex(), cascade.mXi(), cascade.pt(), cascade.yXi(), cascade.phi());
      } else if constexpr (cascadeType == kOmega || cascadeType == kAntiOmega) {
        histos.fill(HIST("h3OmegaCandidatesMass"), collision.centFT0M(), cascade.pt(), cascade.mOmega());
        omegaReducedCandidatesData(collision.globalIndex(), cascade.mOmega(), cascade.pt(), cascade.yOmega(), cascade.phi());
      }
    }
  }

  template <bool isMC, CascadeType... cascadeTypes, typename TColl, typename TCascade>
  void evaluateAllHypotheses(TColl const& collision, TCascade const& cascade, std::optional<std::reference_wrapper<const aod::McParticles>> mcParticlesOpt = std::nullopt)
  {
    (processCascadeType<isMC, cascadeTypes>(collision, cascade, mcParticlesOpt), ...);
  }

  void
    processData(FilteredSelCollisions::iterator const& collision, FullCascades const& cascades, DauTracks const&)
  {
    for (const auto& cascade : cascades) {
      evaluateAllHypotheses<false, kXi, kAntiXi, kOmega, kAntiOmega>(collision, cascade);
    }
  }

  PROCESS_SWITCH(CascadeReducedCandProducer, processData, "Process function to select reduced Xi candidates in Data or in McReco (w/o McTruth) analysis", true);

  void processMCReco(FilteredSimCollisions::iterator const& collision, FullMCCascades const& cascades, DauMCTracks const&, aod::McParticles const& mcParticles)
  {
    for (const auto& cascade : cascades) {
      evaluateAllHypotheses<true, kXi, kAntiXi, kOmega, kAntiOmega>(collision, cascade, mcParticles);
    }
  }

  PROCESS_SWITCH(CascadeReducedCandProducer, processMCReco, "Process function to select reduced Xi candidates in MCReco w MC truth", false);
};

struct PionTrackProducer {
  // Produce the table with the pion tracks information
  Produces<aod::PionTracksData> pionTracksData;
  Produces<aod::PionTracksMcReco> pionTracksMcReco;

  HistogramRegistry histos{"pionTracks", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for analysis mode
  // Configurable<int> analysisMode{"analysisMode", 1, "Analysis mode: 0 - old method with online normalization, 1 - new method with correlations"};

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables for pion track selection
  struct : ConfigurableGroup {
    Configurable<float> cMinPionPtcut{"cMinPionPtcut", 0.2f, "Track minimum pt cut"};
    Configurable<float> nSigmaCutTPCPrimPion{"nSigmaCutTPCPrimPion", 2.0f, "Value of the TPC Nsigma cut for primary Pions"};
    Configurable<float> nSigmaCutCombinedPi{"nSigmaCutCombinedPi", 2.0f, "Value of the TPC and TOF Nsigma cut for Pions"};

    Configurable<bool> cfgIsTOFChecked{"cfgIsTOFChecked", true, "Is TOF checked in PID for pions"};
    Configurable<std::vector<float>> cMaxDCArToPVPion{"cMaxDCArToPVPion", {0.004f, 0.013f, 1.0f}, "Track DCAr cut to PV for Pions"};
    Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0f, "Track DCAz cut to PV Maximum"};
    Configurable<bool> cfgIsDCAzParameterized{"cfgIsDCAzParameterized", false, "IsDCAzParameterized"};
    Configurable<std::vector<float>> cMaxDCAzToPVPion{"cMaxDCAzToPVPion", {0.004f, 0.013f, 1.0f}, "Track DCAz cut to PV for Pions"};

    Configurable<int> minTPCnClsFound{"minTPCnClsFound", 70, "min number of found TPC clusters"};
    Configurable<int> minITSnCls{"minITSnCls", 4, "min number of ITS clusters"};

    // Configurable<bool> forceTOF{"forceTOF", false, "force the TOF signal for the PID"};
    Configurable<float> tofPIDThreshold{"tofPIDThreshold", 0.5, "minimum pT after which TOF PID is applicable"};
    Configurable<std::vector<int>> trkPIDspecies{"trkPIDspecies", std::vector<int>{o2::track::PID::Kaon, o2::track::PID::Proton}, "Trk sel: Particles species for PID rejection, kaon, proton"};
    Configurable<bool> applyElRejection{"applyElRejection", false, "Apply or not the electron rejection"};
    Configurable<std::pair<float, float>> pidRangeTPCEl{"pidRangeTPCEl", {-3.0f, 5.0f}, "nSigma TPC range for electrons"};
    Configurable<float> pidTPCMaxHadrons{"pidTPCMaxHadrons", 3.0f, "maximum nSigma TPC for hadrons"};
    Configurable<float> pidTOFMaxHadrons{"pidTOFMaxHadrons", 3.0f, "maximum nSigma TOF for hadrons"};
    Configurable<std::vector<float>> pidTPCMax{"pidTPCMax", std::vector<float>{3.0f, 3.0f}, "maximum nSigma TPC"};
    Configurable<std::vector<float>> pidTOFMax{"pidTOFMax", std::vector<float>{3.0f, 3.0f}, "maximum nSigma TOF"};

    Configurable<float> cfgYAcceptance{"cfgYAcceptance", 0.5f, "Rapidity acceptance"};
  } trackConfigs;

  // Configurable on pion pT bins
  Configurable<std::vector<double>> binspTPi{"binspTPi", {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0}, "pT bin limits for pions"};

  static constexpr double MassPi = o2::constants::physics::MassPiPlus;

  // Filter on default selected collisions
  // Filter collisionFilter = aod::lf_selection_default_collision::defaultSel == true;
  Filter collisionFilter = aod::lf_selection_event::defaultSel == true;

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::PVMults, aod::PhiStrangeEvtSelDataLike>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;

  using FilteredSelCollisions = soa::Filtered<SelCollisions>;
  using FilteredSimCollisions = soa::Filtered<SimCollisions>;

  // Defining the type of the tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;

  void init(InitContext&)
  {
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};
    AxisSpec binnedpTPiAxis{(std::vector<double>)binspTPi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec nSigmaPiAxis = {100, -10.0f, 10.0f, "N#sigma #pi"};

    histos.add("h3PionTPCnSigma", "Pion TPC nSigma distribution", kTH3F, {binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis});
    histos.add("h3PionTOFnSigma", "Pion TOF nSigma distribution", kTH3F, {binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis});

    histos.add("h2RecMCDCAxyPrimPi", "Dcaxy distribution vs pt for Primary Pions", kTH2F, {binnedpTPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    histos.add("h2RecMCDCAxySecWeakDecayPi", "Dcaz distribution vs pt for Secondary Pions from Weak Decay", kTH2F, {binnedpTPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    histos.add("h2RecMCDCAxySecMaterialPi", "Dcaxy distribution vs pt for Secondary Pions from Material", kTH2F, {binnedpTPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
  }

  /*// PID selection for Pions
  template <typename T>
  bool pidSelectionPion(const T& track)
  {
    for (size_t speciesIndex = 0; speciesIndex < trackConfigs.trkPIDspecies->size(); ++speciesIndex) {
      auto const& pid = trackConfigs.trkPIDspecies->at(speciesIndex);
      auto nSigmaTPC = aod::pidutils::tpcNSigma(pid, track);

      if (trackConfigs.forceTOF && !track.hasTOF()) {
        return false;
      }

      if (speciesIndex == 0) { // First species logic
        if (std::abs(nSigmaTPC) >= trackConfigs.pidTPCMax->at(speciesIndex)) {
          return false; // TPC check failed
        }
        if (trackConfigs.forceTOF || (track.pt() >= trackConfigs.tofPIDThreshold && track.hasTOF())) {
          auto nSigmaTOF = aod::pidutils::tofNSigma(pid, track);
          if (std::abs(nSigmaTOF) >= trackConfigs.pidTOFMax->at(speciesIndex)) {
            return false; // TOF check failed
          }
        }
      } else {                                                                // Other species logic
        if (std::abs(nSigmaTPC) < trackConfigs.pidTPCMax->at(speciesIndex)) { // Check TPC nSigma  first
          if (track.hasTOF()) {
            auto nSigmaTOF = aod::pidutils::tofNSigma(pid, track);
            if (std::abs(nSigmaTOF) < trackConfigs.pidTOFMax->at(speciesIndex)) {
              return false; // Reject if both TPC and TOF are within thresholds
            }
          } else {
            return false; // Reject if only TPC is within threshold and TOF is unavailable
          }
        }
      }
    }

    return true;
  }*/

  // PID Hypotheses rejection for Pions
  template <typename T>
  bool pidHypothesesRejection(const T& track)
  {
    // Electron rejection (if enabled)
    if (trackConfigs.applyElRejection) {
      auto nSigmaTPCEl = aod::pidutils::tpcNSigma(o2::track::PID::Electron, track);

      if (nSigmaTPCEl > trackConfigs.pidRangeTPCEl->first && nSigmaTPCEl < trackConfigs.pidRangeTPCEl->second) {
        auto nSigmaTPCPi = aod::pidutils::tpcNSigma(o2::track::PID::Pion, track);
        auto nSigmaTPCKa = aod::pidutils::tpcNSigma(o2::track::PID::Kaon, track);
        auto nSigmaTPCPr = aod::pidutils::tpcNSigma(o2::track::PID::Proton, track);

        if (std::abs(nSigmaTPCPi) > trackConfigs.pidTPCMaxHadrons &&
            std::abs(nSigmaTPCKa) > trackConfigs.pidTPCMaxHadrons &&
            std::abs(nSigmaTPCPr) > trackConfigs.pidTPCMaxHadrons) {
          return false;
        }
      }
    }

    // Other hadron species rejection
    for (const auto& pid : {o2::track::PID::Kaon, o2::track::PID::Proton}) {
      auto nSigmaTPC = aod::pidutils::tpcNSigma(pid, track);

      if (std::abs(nSigmaTPC) < trackConfigs.pidTPCMaxHadrons) {
        // Reject if only TPC is within threshold and TOF is unavailable or if both TPC and TOF are within thresholds
        if (!track.hasTOF() || std::abs(aod::pidutils::tofNSigma(pid, track)) < trackConfigs.pidTOFMaxHadrons) {
          return false;
        }
      }
    }

    // Other hadron species rejection
    /*for (size_t speciesIndex = 0; speciesIndex < trackConfigs.trkPIDspecies->size(); ++speciesIndex) {
      auto const& pid = trackConfigs.trkPIDspecies->at(speciesIndex);
      auto nSigmaTPC = aod::pidutils::tpcNSigma(pid, track);

      if (std::abs(nSigmaTPC) < trackConfigs.pidTPCMax->at(speciesIndex)) { // Check TPC nSigma first
        if (track.hasTOF()) {
          auto nSigmaTOF = aod::pidutils::tofNSigma(pid, track);
          if (std::abs(nSigmaTOF) < trackConfigs.pidTOFMax->at(speciesIndex)) {
            return false; // Reject if both TPC and TOF are within thresholds
          }
        } else {
          return false; // Reject if only TPC is within threshold and TOF is unavailable
        }
      }
    }*/

    return true;
  }

  // Track selection for Pions
  template <typename T>
  bool selectionPion(const T& track)
  {
    if (!track.isGlobalTrackWoDCA()) {
      return false;
    }

    if (track.itsNCls() < trackConfigs.minITSnCls) {
      return false;
    }
    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound) {
      return false;
    }

    if (track.pt() < trackConfigs.cMinPionPtcut) {
      return false;
    }

    if (std::abs(track.dcaXY()) > trackConfigs.cMaxDCArToPVPion->at(0) + (trackConfigs.cMaxDCArToPVPion->at(1) / std::pow(track.pt(), trackConfigs.cMaxDCArToPVPion->at(2)))) {
      return false;
    }
    if (trackConfigs.cfgIsDCAzParameterized) {
      if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVPion->at(0) + (trackConfigs.cMaxDCAzToPVPion->at(1) / std::pow(track.pt(), trackConfigs.cMaxDCAzToPVPion->at(2)))) {
        return false;
      }
    } else {
      if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVcut) {
        return false;
      }
    }

    if (trackConfigs.cfgIsTOFChecked && track.pt() >= trackConfigs.tofPIDThreshold && !track.hasTOF()) {
      return false;
    }

    /*if (analysisMode == 1 && !pidSelectionPion(track))
      return false;*/

    if (!pidHypothesesRejection(track)) {
      return false;
    }

    /*
    if (analysisMode == 1) {
      if (track.pt() < trackConfigs.tofPIDThreshold && std::abs(track.tpcNSigmaPi()) >= trackConfigs.nSigmaCutTPCPrimPion)
        return false;
      if (trackConfigs.cfgIsTOFChecked && track.pt() >= trackConfigs.tofPIDThreshold && (std::pow(track.tofNSigmaPi(), 2) + std::pow(track.tpcNSigmaPi(), 2)) >= std::pow(trackConfigs.nSigmaCutCombinedPi, 2))
        return false;
    }
    */

    if (std::abs(track.rapidity(MassPi)) > trackConfigs.cfgYAcceptance) {
      return false;
    }

    return true;
  }

  // void processData(SelCollisions::iterator const& collision, FullTracks const& fullTracks)
  void processData(FilteredSelCollisions::iterator const& collision, FullTracks const& fullTracks)
  {
    for (const auto& track : fullTracks) {
      if (!selectionPion(track)) {
        continue;
      }

      histos.fill(HIST("h3PionTPCnSigma"), collision.centFT0M(), track.pt(), track.tpcNSigmaPi());
      histos.fill(HIST("h3PionTOFnSigma"), collision.centFT0M(), track.pt(), track.tofNSigmaPi());

      pionTracksData(collision.globalIndex(), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt(), track.rapidity(MassPi), track.phi(), track.hasTOF());
    }
  }

  PROCESS_SWITCH(PionTrackProducer, processData, "Process function to select pion tracks in Data or in McReco (w/o MC truth) analysis", true);

  // void processMCReco(SimCollisions::iterator const& collision, FullMCTracks const& fullTracks, aod::McParticles const& mcParticles)
  void processMCReco(FilteredSimCollisions::iterator const& collision, FullMCTracks const& fullTracks, aod::McParticles const& mcParticles)
  {
    for (const auto& track : fullTracks) {
      if (!selectionPion(track)) {
        continue;
      }
      if (!track.has_mcParticle()) {
        continue;
      }

      const auto trackMcParticle = mcParticles.rawIteratorAt(track.mcParticleId());
      if (std::abs(trackMcParticle.pdgCode()) != PDG_t::kPiPlus) {
        continue;
      }

      if (trackMcParticle.isPhysicalPrimary()) {
        histos.fill(HIST("h2RecMCDCAxyPrimPi"), track.pt(), track.dcaXY());
      } else {
        if (trackMcParticle.getProcess() == TMCProcess::kPDecay) { // Selection of secondary pions from weak decay
          histos.fill(HIST("h2RecMCDCAxySecWeakDecayPi"), track.pt(), track.dcaXY());
        } else { // Selection of secondary pions from material interactions
          histos.fill(HIST("h2RecMCDCAxySecMaterialPi"), track.pt(), track.dcaXY());
        }
        continue;
      }

      pionTracksMcReco(collision.globalIndex(), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt(), track.rapidity(MassPi), track.phi(), track.hasTOF());
    }
  }

  PROCESS_SWITCH(PionTrackProducer, processMCReco, "Process function to select pion tracks in MCReco w MC truth", false);
};

struct EventSelectionProducer {
  // Produce the table with the event selection information
  Produces<aod::PhiStrangeEvtSelDataLike> phiStrangeEvtSelDataLike;
  Produces<aod::PhiStrangeEvtSelMcGen> phiStrangeEvtSelMcGen;

  // Produces<aod::PhiStrangeDefEvtSelDataLike> phiStrangeDefEvtSelDataLike;
  // Produces<aod::PhiStrangeDefEvtSelMcGen> phiStrangeDefEvtSelMcGen;

  /*Produces<aod::DefaultSelectionData> defaultSelectionData;
  Produces<aod::DefaultSelectionMcGen> defaultSelectionMcGen;

  Produces<aod::PhimesonSelectionData> phimesonSelectionData;
  Produces<aod::PhimesonSelectionMcGen> phimesonSelectionMcGen;*/

  HistogramRegistry histos{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for event selection
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurable for MC
  Configurable<bool> cfgiskNoITSROFrameBorder{"cfgiskNoITSROFrameBorder", false, "kNoITSROFrameBorder request on MC collisions"};

  // Configurables on phi selection
  struct : ConfigurableGroup {
    Configurable<float> minMPhiSignal{"minMPhiSignal", 1.0095f, "Upper limits on Phi mass for signal extraction"};
    Configurable<float> maxMPhiSignal{"maxMPhiSignal", 1.029f, "Upper limits on Phi mass for signal extraction"};
  } phiConfigs;

  Configurable<bool> withTrueGenPhi{"withTrueGenPhi", false, "Require the presence of a true generated phi meson in the event"};

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;
  using MCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  // Defining the type of the phi's daughter tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;

  // Necessary service to flag INEL>0 events in GenMC
  Service<o2::framework::O2DatabasePDG> pdgDB{};

  void init(InitContext&)
  {
    // Defining histogram axes
    AxisSpec vertexZAxis = {100, -cutZVertex, cutZVertex, "vrtx_{Z} [cm]"};
    AxisSpec multAxis = {120, 0.0f, 120.0f, "centFT0M"};
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};

    // Booking histograms for event selection QA
    // Number of events per selection in Data
    histos.add("hEventSelectionData", "hEventSelectionData", kTH1F, {{5, -0.5f, 4.5f}});
    histos.get<TH1>(HIST("hEventSelectionData"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelectionData"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos.get<TH1>(HIST("hEventSelectionData"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    histos.get<TH1>(HIST("hEventSelectionData"))->GetXaxis()->SetBinLabel(4, "INEL>0 cut");
    histos.get<TH1>(HIST("hEventSelectionData"))->GetXaxis()->SetBinLabel(5, "With at least a #phi cand");

    // Number of MC events per selection in MC
    histos.add("hEventSelectionMC", "hEventSelectionMC", kTH1F, {{8, -0.5f, 7.5f}});
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(2, "kIsTriggerTVX");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(3, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(5, "posZ cut");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(6, "INEL>0 cut");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(7, "With at least a gen coll");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(8, "With at least a #phi cand");

    // Event information
    histos.add("hVertexZ", "Vertex Z", kTH1F, {vertexZAxis});
    histos.add("hVertexZWPhi", "Vertex Z with a Phi Candidate", kTH1F, {vertexZAxis});
    histos.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {multAxis});
    histos.add("hMultiplicityPercentWPhi", "Multiplicity Percentile in Events with a Phi Candidate", kTH1F, {multAxis});
    histos.add("h2VertexZvsMult", "Vertex Z vs Multiplicity Percentile", kTH2F, {vertexZAxis, binnedmultAxis});
    histos.add("h2VertexZvsMultWPhi", "Vertex Z vs Multiplicity Percentile with a Phi Candidate", kTH2F, {vertexZAxis, binnedmultAxis});
  }

  // Default event selection
  template <bool isMC, typename MC_T = void, typename T>
  bool defaultEventSelection(const T& collision)
  {
    float multPercentile{0.0f};

    if constexpr (!isMC) {                         // data event
      histos.fill(HIST("hEventSelectionData"), 0); // all collisions
      if (!collision.sel8()) {
        return false;
      }
      histos.fill(HIST("hEventSelectionData"), 1); // sel8 collisions
      if (std::abs(collision.posZ()) >= cutZVertex) {
        return false;
      }
      histos.fill(HIST("hEventSelectionData"), 2); // vertex-Z selected
      if (!collision.isInelGt0()) {
        return false;
      }
      histos.fill(HIST("hEventSelectionData"), 3); // INEL>0 collisions

      multPercentile = collision.centFT0M();
    } else { // MCreco event
      static_assert(!std::is_same_v<MC_T, void>, "Need to set MC_T to MCCollisions for isMC = true");

      histos.fill(HIST("hEventSelectionMC"), 0); // all collisions
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        return false;
      }
      histos.fill(HIST("hEventSelectionMC"), 1); // kIsTriggerTVX collisions
      if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        return false;
      }
      histos.fill(HIST("hEventSelectionMC"), 2); // kNoTimeFrameBorder collisions
      if (cfgiskNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        return false;
      }
      histos.fill(HIST("hEventSelectionMC"), 3); // kNoITSROFrameBorder collisions (by default not requested by the selection)
      if (std::abs(collision.posZ()) > cutZVertex) {
        return false;
      }
      histos.fill(HIST("hEventSelectionMC"), 4); // vertex-Z selected
      if (!collision.isInelGt0()) {
        return false;
      }
      histos.fill(HIST("hEventSelectionMC"), 5); // INEL>0 collisions
      if (!collision.has_mcCollision()) {
        return false;
      }
      histos.fill(HIST("hEventSelectionMC"), 6); // with at least a gen collision

      const auto& mcCollision = collision.template mcCollision_as<MC_T>();
      multPercentile = mcCollision.centFT0M();
    }

    histos.fill(HIST("hVertexZ"), collision.posZ());
    histos.fill(HIST("hMultiplicityPercent"), multPercentile);
    histos.fill(HIST("h2VertexZvsMult"), collision.posZ(), multPercentile);

    return true;
  }

  // Check if the event has at least one phi candidate
  template <bool isMC, typename MC_T = void, typename T1, typename T2>
  bool eventHasPhi(const T1& collision, const T2& phiCandidates)
  {
    if (withTrueGenPhi) {
      if (phiCandidates.size() < 1) {
        return false;
      }
    } else {
      uint16_t nPhi{0};

      for (const auto& phiCand : phiCandidates) {
        if (phiCand.inMassRegion(phiConfigs.minMPhiSignal, phiConfigs.maxMPhiSignal)) {
          nPhi++;
        }

        // histos.fill(HIST("hEta"), track1.eta());
        // histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcInnerParam(), track1.tpcNSigmaKa());
        // histos.fill(HIST("hNsigmaKaonTOF"), track1.tpcInnerParam(), track1.tofNSigmaKa());
        // histos.fill(HIST("h2DauTracksPhiDCAxy"), track1.pt(), track1.dcaXY());
        // histos.fill(HIST("h2DauTracksPhiDCAz"), track1.pt(), track1.dcaZ());
      }

      if (nPhi == 0) {
        return false;
      }
    }

    float multPercentile{0.0f};

    if constexpr (!isMC) {
      histos.fill(HIST("hEventSelectionData"), 4);

      multPercentile = collision.centFT0M();
    } else {
      if constexpr (std::is_same_v<MC_T, void>) {
        multPercentile = collision.centFT0M();
      } else {
        histos.fill(HIST("hEventSelectionMC"), 7);

        const auto& mcCollision = collision.template mcCollision_as<MC_T>();
        multPercentile = mcCollision.centFT0M();
      }
    }

    histos.fill(HIST("hVertexZWPhi"), collision.posZ());
    histos.fill(HIST("hMultiplicityPercentWPhi"), multPercentile);
    histos.fill(HIST("h2VertexZvsMultWPhi"), collision.posZ(), multPercentile);

    return true;
  }

  // Default event selection
  void processData(SelCollisions::iterator const& collision)
  {
    // phiStrangeDefEvtSelDataLike(defaultEventSelection<false>(collision));
    phiStrangeEvtSelDataLike(defaultEventSelection<false>(collision), false);
    // defaultSelectionData(defaultEventSelection<false>(collision), false);
  }

  PROCESS_SWITCH(EventSelectionProducer, processData, "Process function to select default events in Data", true);

  void processMCReco(SimCollisions::iterator const& collision, MCCollisions const&)
  {
    // phiStrangeDefEvtSelDataLike(defaultEventSelection<true, MCCollisions>(collision));
    phiStrangeEvtSelDataLike(defaultEventSelection<true, MCCollisions>(collision), false);
    // defaultSelectionData(defaultEventSelection<true, MCCollisions>(collision));
  }

  PROCESS_SWITCH(EventSelectionProducer, processMCReco, "Process function to select default events in MCReco", false);

  void processMCGen(MCCollisions::iterator const&, aod::McParticles const& mcParticles)
  {

    // phiStrangeDefEvtSelMcGen(pwglf::isINELgt0mc(mcParticles, pdgDB));
    phiStrangeEvtSelMcGen(pwglf::isINELgt0mc(mcParticles, pdgDB), false);
    // defaultSelectionMcGen(pwglf::isINELgt0mc(mcParticles, pdgDB));
  }

  PROCESS_SWITCH(EventSelectionProducer, processMCGen, "Process function to select default events in MCGen", false);

  // Default event selection + phi meson requirement
  void processDataWPhi(SelCollisions::iterator const& collision, aod::PhimesonCandidatesData const& phiCandidatesData)
  {
    phiStrangeEvtSelDataLike(defaultEventSelection<false>(collision), eventHasPhi<false>(collision, phiCandidatesData));
    // phimesonSelectionData(defaultEventSelection<false>(collision) && eventHasPhi<false>(collision, phiCandidatesData));
  }

  PROCESS_SWITCH(EventSelectionProducer, processDataWPhi, "Process function to select events with Phi mesons in Data", false);

  void processMCRecoWPhi(SimCollisions::iterator const& collision, MCCollisions const&, aod::PhimesonCandidatesData const& phiCandidatesData)
  {
    phiStrangeEvtSelDataLike(defaultEventSelection<true, MCCollisions>(collision), eventHasPhi<true, MCCollisions>(collision, phiCandidatesData));
    // phimesonSelectionData(defaultEventSelection<true, MCCollisions>(collision) && eventHasPhi<true, MCCollisions>(collision, phiCandidatesData));
  }

  PROCESS_SWITCH(EventSelectionProducer, processMCRecoWPhi, "Process function to select events with Phi mesons in MCReco", false);

  void processMCGenWPhi(MCCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles, aod::PhimesonCandidatesMcGen const& phiCandidatesMcGen)
  {
    phiStrangeEvtSelMcGen(pwglf::isINELgt0mc(mcParticles, pdgDB), eventHasPhi<true>(mcCollision, phiCandidatesMcGen));
    // phimesonSelectionMcGen(pwglf::isINELgt0mc(mcParticles, pdgDB) && eventHasPhi<true>(mcCollision, phiCandidatesMcGen));
  }

  PROCESS_SWITCH(EventSelectionProducer, processMCGenWPhi, "Process function to select events with Phi mesons in MCGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhiMesonCandProducer>(cfgc),
                      adaptAnalysisTask<V0ReducedCandProducer>(cfgc),
                      adaptAnalysisTask<CascadeReducedCandProducer>(cfgc),
                      adaptAnalysisTask<PionTrackProducer>(cfgc),
                      adaptAnalysisTask<EventSelectionProducer>(cfgc)};
}
