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

#include "PWGLF/DataModel/LFPhiStrangeCorrelationTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/TableHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <TMath.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
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

  HistogramRegistry histos{"phiCandidates", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

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
    Configurable<float> maxMPhi{"maxMPhi", 1.5f, "Maximum mass for Phi candidates"};
    Configurable<float> minPhiPt{"minPhiPt", 0.4f, "Minimum pT for Phi candidates"};
    Configurable<float> cfgYAcceptance{"cfgYAcceptance", 0.5f, "Rapidity acceptance"};
  } phiConfigs;

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;

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
    AxisSpec massPhiAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};

    histos.add("h1PhiCandidateMass", "Phi candidate invariant mass", kTH1F, {massPhiAxis});
  }

  // Topological track selection
  template <typename T>
  bool selectionTrackResonance(const T& track)
  {
    if (trackConfigs.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (trackConfigs.cfgPVContributor && !track.isPVContributor())
      return false;

    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;

    if (track.pt() < trackConfigs.cMinKaonPtcut)
      return false;

    if (std::abs(track.dcaXY()) > trackConfigs.cMaxDCArToPVPhi->at(0) + (trackConfigs.cMaxDCArToPVPhi->at(1) / std::pow(track.pt(), trackConfigs.cMaxDCArToPVPhi->at(2))))
      return false;
    if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVcut)
      return false;

    return true;
  }

  // PIDKaon track selection
  template <typename T>
  bool selectionPIDKaonpTdependent(const T& track)
  {
    if (track.pt() < trackConfigs.pTToUseTOF && std::abs(track.tpcNSigmaKa()) < trackConfigs.nSigmaCutTPCKa)
      return true;
    if (track.pt() >= trackConfigs.pTToUseTOF && track.hasTOF() && (std::pow(track.tofNSigmaKa(), 2) + std::pow(track.tpcNSigmaKa(), 2)) < std::pow(trackConfigs.nSigmaCutCombinedKa, 2))
      return true;

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

  void processData(SelCollisions::iterator const& collision, FullTracks const&)
  {
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
        continue;

      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
          continue;

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);

        if (recPhi.Pt() < phiConfigs.minPhiPt)
          continue;
        if (recPhi.M() > phiConfigs.maxMPhi)
          continue;
        if (std::abs(recPhi.Rapidity()) > phiConfigs.cfgYAcceptance)
          continue;

        histos.fill(HIST("h1PhiCandidateMass"), recPhi.M());

        phimesonCandidatesData(collision.globalIndex(), recPhi.M(), recPhi.Pt(), recPhi.Rapidity(), recPhi.Phi());
      }
    }
  }

  PROCESS_SWITCH(PhiMesonCandProducer, processData, "Process function to select Phi meson candidates in Data or in ", true);

  void processMCRecoDataLike(SimCollisions::iterator const& collision, FullMCTracks const&)
  {
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
        continue;

      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
          continue;

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);

        if (recPhi.Pt() < phiConfigs.minPhiPt)
          continue;
        if (recPhi.M() > phiConfigs.maxMPhi)
          continue;
        if (std::abs(recPhi.Rapidity()) > phiConfigs.cfgYAcceptance)
          continue;

        phimesonCandidatesData(collision.globalIndex(), recPhi.M(), recPhi.Pt(), recPhi.Rapidity(), recPhi.Phi());
      }
    }
  }

  PROCESS_SWITCH(PhiMesonCandProducer, processMCRecoDataLike, "Process function to select Phi meson candidates in MCReco w/o MC truth", false);

  void processMCReco(SimCollisions::iterator const& collision, FullMCTracks const&, aod::McParticles const& mcParticles)
  {
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
        continue;
      if (!track1.has_mcParticle())
        continue;
      const auto track1McParticle = mcParticles.rawIteratorAt(track1.mcParticleId());
      if (track1McParticle.pdgCode() != PDG_t::kKPlus || !track1McParticle.isPhysicalPrimary())
        continue;

      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
          continue;
        if (!track2.has_mcParticle())
          continue;
        const auto track2McParticle = mcParticles.rawIteratorAt(track2.mcParticleId());
        if (track2McParticle.pdgCode() != PDG_t::kKMinus || !track2McParticle.isPhysicalPrimary())
          continue;

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);

        if (recPhi.Pt() < phiConfigs.minPhiPt)
          continue;
        if (recPhi.M() > phiConfigs.maxMPhi)
          continue;
        if (std::abs(recPhi.Rapidity()) > phiConfigs.cfgYAcceptance)
          continue;

        const auto track1mcPartMotherIndexes = track1McParticle.mothersIds();
        const auto track2mcPartMotherIndexes = track2McParticle.mothersIds();

        auto genPhiMaybe = [&]() -> std::optional<aod::McParticles::iterator> {
          for (const auto& mother1Index : track1mcPartMotherIndexes) {
            for (const auto& mother2Index : track2mcPartMotherIndexes) {
              if (mother1Index != mother2Index)
                continue;

              const auto motherMcParticle = mcParticles.rawIteratorAt(mother1Index);
              if (std::abs(motherMcParticle.pdgCode()) == o2::constants::physics::Pdg::kPhi)
                return motherMcParticle;
            }
          }

          return std::nullopt;
        }();

        if (!genPhiMaybe)
          continue;
        const auto genPhi = *genPhiMaybe;

        phimesonCandidatesMcReco(collision.globalIndex(), recPhi.M(), genPhi.pt(), genPhi.y(), genPhi.phi());
      }
    }
  }

  PROCESS_SWITCH(PhiMesonCandProducer, processMCReco, "Process function to select Phi meson candidates in MCReco w MC truth", false);

  void processMCGen(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    for (const auto& mcParticle1 : mcParticles) {
      if (!mcParticle1.isPhysicalPrimary() || std::abs(mcParticle1.eta()) > trackConfigs.etaMax)
        continue;

      for (const auto& mcParticle2 : mcParticles) {
        if (!mcParticle2.isPhysicalPrimary() || std::abs(mcParticle2.eta()) > trackConfigs.etaMax)
          continue;

        if (!(mcParticle1.pdgCode() == PDG_t::kKPlus && mcParticle2.pdgCode() == PDG_t::kKMinus) &&
            !(mcParticle1.pdgCode() == PDG_t::kKMinus && mcParticle2.pdgCode() == PDG_t::kKPlus))
          continue;

        ROOT::Math::PxPyPzMVector genKPair = recMother(mcParticle1, mcParticle2, massKa, massKa);

        if (genKPair.Pt() < phiConfigs.minPhiPt)
          continue;
        if (genKPair.M() > phiConfigs.maxMPhi)
          continue;
        if (std::abs(genKPair.Rapidity()) > phiConfigs.cfgYAcceptance)
          continue;

        phimesonCandidatesMcGen(mcCollision.globalIndex(), genKPair.M(), genKPair.Pt(), genKPair.Rapidity(), genKPair.Phi());
      }
    }
  }

  PROCESS_SWITCH(PhiMesonCandProducer, processMCGen, "Process function to select Phi meson candidates in MCGen", false);
};

struct PhiMesonSelCollision {
  // Produce the table with the event selection information
  Produces<aod::PhimesonSelectionData> phimesonSelectionData;
  Produces<aod::PhimesonSelectionMcGen> phimesonSelectionMcGen;

  HistogramRegistry histos{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for selection type
  Configurable<int> selectionType{"selectionType", 1, "Selection type: 0 - default selection only, 1 - default + phi meson selection"};

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

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;
  using MCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  // Defining the type of the phi's daughter tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;

  // Necessary service to flag INEL>0 events in GenMC
  Service<o2::framework::O2DatabasePDG> pdgDB;

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

    // Phi's daughter tracks information
    /*histos.add("hEta", "Eta of Kaon candidates", kTH1F, {{100, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});
    histos.add("h2DauTracksPhiDCAxy", "DCAxy distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    histos.add("h2DauTracksPhiDCAz", "DCAz distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});*/
  }

  // Default event selection
  template <bool isMC, typename MC_T = void, typename T>
  bool defaultEventSelection(const T& collision)
  {
    float multPercentile{0.0f};

    if constexpr (!isMC) {                         // data event
      histos.fill(HIST("hEventSelectionData"), 0); // all collisions
      if (!collision.sel8())
        return false;
      histos.fill(HIST("hEventSelectionData"), 1); // sel8 collisions
      if (std::abs(collision.posZ()) >= cutZVertex)
        return false;
      histos.fill(HIST("hEventSelectionData"), 2); // vertex-Z selected
      if (!collision.isInelGt0())
        return false;
      histos.fill(HIST("hEventSelectionData"), 3); // INEL>0 collisions

      multPercentile = collision.centFT0M();
    } else { // MCreco event
      static_assert(!std::is_same_v<MC_T, void>, "Need to set MC_T to MCCollisions for isMC = true");

      histos.fill(HIST("hEventSelectionMC"), 0); // all collisions
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX))
        return false;
      histos.fill(HIST("hEventSelectionMC"), 1); // kIsTriggerTVX collisions
      if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        return false;
      histos.fill(HIST("hEventSelectionMC"), 2); // kNoTimeFrameBorder collisions
      if (cfgiskNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
        return false;
      histos.fill(HIST("hEventSelectionMC"), 3); // kNoITSROFrameBorder collisions (by default not requested by the selection)
      if (std::abs(collision.posZ()) > cutZVertex)
        return false;
      histos.fill(HIST("hEventSelectionMC"), 4); // vertex-Z selected
      if (!collision.isInelGt0())
        return false;
      histos.fill(HIST("hEventSelectionMC"), 5); // INEL>0 collisions
      if (!collision.has_mcCollision())
        return false;
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
    uint16_t nPhi{0};

    for (const auto& phiCand : phiCandidates) {

      if (phiCand.inMassRegion(phiConfigs.minMPhiSignal, phiConfigs.maxMPhiSignal))
        nPhi++;

      // histos.fill(HIST("hEta"), track1.eta());
      // histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcInnerParam(), track1.tpcNSigmaKa());
      // histos.fill(HIST("hNsigmaKaonTOF"), track1.tpcInnerParam(), track1.tofNSigmaKa());
      // histos.fill(HIST("h2DauTracksPhiDCAxy"), track1.pt(), track1.dcaXY());
      // histos.fill(HIST("h2DauTracksPhiDCAz"), track1.pt(), track1.dcaZ());
    }

    if (nPhi == 0)
      return false;

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

  void processData(SelCollisions::iterator const& collision, aod::PhimesonCandidatesData const& phiCandidatesData)
  {
    phimesonSelectionData(defaultEventSelection<false>(collision) && selectionType == 1 ? eventHasPhi<false>(collision, phiCandidatesData) : true);
  }

  PROCESS_SWITCH(PhiMesonSelCollision, processData, "Process function to select events with Phi mesons in Data", true);

  void processMCReco(SimCollisions::iterator const& collision, MCCollisions const&, aod::PhimesonCandidatesData const& phiCandidatesData)
  {
    phimesonSelectionData(defaultEventSelection<true, MCCollisions>(collision) && selectionType == 1 ? eventHasPhi<true, MCCollisions>(collision, phiCandidatesData) : true);
  }

  PROCESS_SWITCH(PhiMesonSelCollision, processMCReco, "Process function to select events with Phi mesons in MCReco", false);

  void processMCGen(MCCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles, aod::PhimesonCandidatesMcGen const& phiCandidatesMcGen)
  {
    phimesonSelectionMcGen(pwglf::isINELgt0mc(mcParticles, pdgDB) && selectionType == 1 ? eventHasPhi<true>(mcCollision, phiCandidatesMcGen) : true);
  }

  PROCESS_SWITCH(PhiMesonSelCollision, processMCGen, "Process function to select events with Phi mesons in MCGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhiMesonCandProducer>(cfgc),
                      adaptAnalysisTask<PhiMesonSelCollision>(cfgc)};
}
