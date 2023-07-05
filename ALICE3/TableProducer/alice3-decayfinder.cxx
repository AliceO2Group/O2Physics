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
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//   Decay finder task for ALICE 3
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Uses specific ALICE 3 PID and performance for studying
//    HF decays. Work in progress: use at your own risk!
//

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/RICH.h"
#include "ALICE3/DataModel/A3DecayFinderTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// simple checkers
#define biton(var, nbit) ((var) |= (static_cast<uint32_t>(1) << (nbit)))
#define bitoff(var, nbit) ((var) &= ~(static_cast<uint32_t>(1) << (nbit))) //((a) &= ~(1ULL<<(b)))
#define bitcheck(var, nbit) ((var) & (static_cast<uint32_t>(1) << (nbit)))

using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov>;

// For MC association in pre-selection
using labeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using tofTracks = soa::Join<aod::Tracks, aod::UpgradeTofs>;
using richTracks = soa::Join<aod::Tracks, aod::RICHs>;
using alice3tracks = soa::Join<aod::Tracks, aod::TracksCov, aod::Alice3DecayMaps, aod::McTrackLabels>;

struct alice3decayPreselector {
  Produces<aod::Alice3DecayMaps> a3decayMaps;

  // Operation and minimisation criteria
  Configurable<float> nSigmaTOF{"nSigmaTOF", 4.0f, "Nsigma for TOF PID (if enabled)"};
  Configurable<float> nSigmaRICH{"nSigmaRICH", 4.0f, "Nsigma for RICH PID (if enabled)"};

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> fitter;

  // for bit-packed maps
  std::vector<uint32_t> selectionMap;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  /// function to check PDG + PDG mother
  template <typename TTrack>
  bool checkPDG(TTrack const& track, int pdgMother, int pdg)
  {
    bool returnValue = false;
    // Association check
    if (track.has_mcParticle()) {
      auto mcParticle = track.template mcParticle_as<aod::McParticles>();
      if (mcParticle.has_mothers()) {
        for (auto& mcParticleMother : mcParticle.template mothers_as<aod::McParticles>()) {
          if (mcParticle.pdgCode() == pdg && mcParticleMother.pdgCode() == pdgMother)
            returnValue = true;
        }
      }
    } // end association check
    return returnValue;
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  void init(InitContext& context)
  {
    // future dev if needed
    histos.add("hToggle", "hToggle", kTH1F, {{10, -0.5f, 9.5f}});
  }

  // go declarative: use partitions instead of "if", then just toggle bits to allow for mask selection later
  Partition<tofTracks> pInnerTOFPi = nabs(aod::upgrade_tof::nSigmaPionInnerTOF) > nSigmaTOF;
  Partition<tofTracks> pInnerTOFKa = nabs(aod::upgrade_tof::nSigmaKaonInnerTOF) > nSigmaTOF;
  Partition<tofTracks> pInnerTOFPr = nabs(aod::upgrade_tof::nSigmaProtonInnerTOF) > nSigmaTOF;
  Partition<tofTracks> pOuterTOFPi = nabs(aod::upgrade_tof::nSigmaPionOuterTOF) > nSigmaTOF;
  Partition<tofTracks> pOuterTOFKa = nabs(aod::upgrade_tof::nSigmaKaonOuterTOF) > nSigmaTOF;
  Partition<tofTracks> pOuterTOFPr = nabs(aod::upgrade_tof::nSigmaProtonOuterTOF) > nSigmaTOF;
  Partition<richTracks> pRICHPi = nabs(aod::alice3rich::richNsigmaPi) > nSigmaRICH;
  Partition<richTracks> pRICHKa = nabs(aod::alice3rich::richNsigmaKa) > nSigmaRICH;
  Partition<richTracks> pRICHPr = nabs(aod::alice3rich::richNsigmaPr) > nSigmaRICH;

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  /// Initialization of mask vectors if uninitialized
  void initializeMasks(int size)
  {
    selectionMap.clear();
    selectionMap.resize(size, 0xFFFFFFFF); // all bits 1, please
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  /// This process function ensures that all V0s are built. It will simply tag everything as true.
  void processInitialize(aod::Tracks const& tracks)
  {
    initializeMasks(tracks.size());
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFilterInnerTOF(tofTracks const& tracks)
  {
    for (auto const& track : pInnerTOFPi)
      bitoff(selectionMap[track.globalIndex()], kInnerTOFPion);
    for (auto const& track : pInnerTOFKa)
      bitoff(selectionMap[track.globalIndex()], kInnerTOFKaon);
    for (auto const& track : pInnerTOFPr)
      bitoff(selectionMap[track.globalIndex()], kInnerTOFProton);
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFilterOuterTOF(tofTracks const& tracks)
  {
    for (auto const& track : pOuterTOFPi)
      bitoff(selectionMap[track.globalIndex()], kOuterTOFPion);
    for (auto const& track : pOuterTOFKa)
      bitoff(selectionMap[track.globalIndex()], kOuterTOFKaon);
    for (auto const& track : pOuterTOFPr)
      bitoff(selectionMap[track.globalIndex()], kOuterTOFProton);
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFilterRICH(richTracks const& tracks)
  {
    for (auto const& track : pRICHPi)
      bitoff(selectionMap[track.globalIndex()], kRICHPion);
    for (auto const& track : pRICHKa)
      bitoff(selectionMap[track.globalIndex()], kRICHKaon);
    for (auto const& track : pRICHPr)
      bitoff(selectionMap[track.globalIndex()], kRICHProton);
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFilterOnMonteCarloTruth(labeledTracks const& tracks, aod::McParticles const&)
  {
    for (auto const& track : tracks) {
      if (!checkPDG(track, 421, -321)) //+421 -> -321 +211
      {
        bitoff(selectionMap[track.globalIndex()], kTrueKaMinusFromD);
        histos.fill(HIST("hToggle"), 0.0);
      } else {
        histos.fill(HIST("hToggle"), 1.0);
      }
      if (!checkPDG(track, -421, +321)) //-421 -> +321 -211
        bitoff(selectionMap[track.globalIndex()], kTrueKaPlusFromD);
      if (!checkPDG(track, 421, +211)) //+421 -> -321 +211
      {
        bitoff(selectionMap[track.globalIndex()], kTruePiPlusFromD);
        histos.fill(HIST("hToggle"), 2.0);
      } else {
        histos.fill(HIST("hToggle"), 3.0);
      }
      if (!checkPDG(track, -421, -211)) //-421 -> +321 -211
        bitoff(selectionMap[track.globalIndex()], kTruePiMinusFromD);
    }
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processPublishDecision(aod::Tracks const& tracks)
  {
    for (uint32_t i = 0; i < tracks.size(); i++) {
      a3decayMaps(selectionMap[i]);
    }
    selectionMap.clear();
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
  PROCESS_SWITCH(alice3decayPreselector, processInitialize, "Initialize (MUST be on)", true);
  PROCESS_SWITCH(alice3decayPreselector, processFilterInnerTOF, "Switch to use inner TOF PID", false);
  PROCESS_SWITCH(alice3decayPreselector, processFilterOuterTOF, "Switch to use outer TOF PID", false);
  PROCESS_SWITCH(alice3decayPreselector, processFilterRICH, "Switch to use RICH", false);
  PROCESS_SWITCH(alice3decayPreselector, processFilterOnMonteCarloTruth, "Switch to use MC truth", false);
  PROCESS_SWITCH(alice3decayPreselector, processPublishDecision, "Fill decision mask table (MUST be on)", true);
  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
};

struct alice3decayFinder {
  SliceCache cache;

  // Operation and minimisation criteria
  Configurable<float> magneticField{"magneticField", 20.0f, "Magnetic field (in kilogauss)"};
  Configurable<bool> doDCAplots{"doDCAplots", true, "do daughter prong DCA plots"};
  Configurable<bool> mcSameMotherCheck{"mcSameMotherCheck", true, "check if tracks come from the same MC mother"};
  Configurable<float> dcaDaughtersSelection{"dcaDaughtersSelection", 1000.0f, "DCA between daughters (cm)"};

  ConfigurableAxis axisEta{"axisEta", {8, -4.0f, +4.0f}, "#eta"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisDCA{"axisDCA", {200, -100, 100}, "DCA (#mum)"};
  ConfigurableAxis axisDMass{"axisDMass", {200, 1.765f, 1.965f}, "D Inv Mass (GeV/c^{2})"};

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> fitter;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Partition<aod::McParticles> trueD = aod::mcparticle::pdgCode == 421;
  Partition<aod::McParticles> trueDbar = aod::mcparticle::pdgCode == -421;

  static constexpr uint32_t trackSelectionPiPlusFromD = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiPlusFromD;
  static constexpr uint32_t trackSelectionPiMinusFromD = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiMinusFromD;
  static constexpr uint32_t trackSelectionKaPlusFromD = 1 << kInnerTOFKaon | 1 << kOuterTOFKaon | 1 << kRICHKaon | 1 << kTrueKaPlusFromD;
  static constexpr uint32_t trackSelectionKaMinusFromD = 1 << kInnerTOFKaon | 1 << kOuterTOFKaon | 1 << kRICHKaon | 1 << kTrueKaMinusFromD;

  Partition<alice3tracks> tracksPiPlusFromD =
    ((aod::a3DecayMap::decayMap & trackSelectionPiPlusFromD) == trackSelectionPiPlusFromD) && aod::track::signed1Pt > 0.0f;

  Partition<alice3tracks> tracksPiMinusFromD =
    ((aod::a3DecayMap::decayMap & trackSelectionPiMinusFromD) == trackSelectionPiMinusFromD) && aod::track::signed1Pt < 0.0f;

  Partition<alice3tracks> tracksKaPlusFromD =
    ((aod::a3DecayMap::decayMap & trackSelectionKaPlusFromD) == trackSelectionKaPlusFromD) && aod::track::signed1Pt > 0.0f;

  Partition<alice3tracks> tracksKaMinusFromD =
    ((aod::a3DecayMap::decayMap & trackSelectionKaMinusFromD) == trackSelectionKaMinusFromD) && aod::track::signed1Pt < 0.0f;

  // Helper struct to pass candidate information
  struct {
    float mass;
    float pt;
    float eta;
  } dmeson;

  template <typename TTrackType>
  bool buildDecayCandidate(TTrackType const& posTrackRow, TTrackType const& negTrackRow, float posMass, float negMass)
  {
    o2::track::TrackParCov posTrack = getTrackParCov(posTrackRow);
    o2::track::TrackParCov negTrack = getTrackParCov(negTrackRow);

    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter.process(posTrack, negTrack);
    } catch (...) {
      return false;
    }
    if (nCand == 0) {
      return false;
    }
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}

    posTrack = fitter.getTrack(0);
    negTrack = fitter.getTrack(1);
    std::array<float, 3> posP;
    std::array<float, 3> negP;
    posTrack.getPxPyPzGlo(posP);
    negTrack.getPxPyPzGlo(negP);

    float dcaDau = TMath::Sqrt(fitter.getChi2AtPCACandidate());
    if (dcaDau > dcaDaughtersSelection)
      return false;

    // return mass
    dmeson.mass = RecoDecay::m(array{array{posP[0], posP[1], posP[2]}, array{negP[0], negP[1], negP[2]}}, array{posMass, negMass});
    dmeson.pt = std::hypot(posP[0] + negP[0], posP[1] + negP[1]);
    dmeson.eta = RecoDecay::eta(array{posP[0] + negP[0], posP[1] + negP[1], posP[2] + negP[2]});
    return true;
  }

  /// function to check if tracks have the same mother in MC
  template <typename TTrackType>
  bool checkSameMother(TTrackType const& track1, TTrackType const& track2)
  {
    bool returnValue = false;
    // Association check
    // There might be smarter ways of doing this in the future
    if (track1.has_mcParticle() && track2.has_mcParticle()) {
      auto mcParticle1 = track1.template mcParticle_as<aod::McParticles>();
      auto mcParticle2 = track2.template mcParticle_as<aod::McParticles>();
      if (mcParticle1.has_mothers() && mcParticle2.has_mothers()) {
        for (auto& mcParticleMother1 : mcParticle1.template mothers_as<aod::McParticles>()) {
          for (auto& mcParticleMother2 : mcParticle2.template mothers_as<aod::McParticles>()) {
            if (mcParticleMother1.globalIndex() == mcParticleMother2.globalIndex()) {
              returnValue = true;
            }
          }
        }
      }
    } // end association check
    return returnValue;
  }

  void init(InitContext& context)
  {
    // initialize O2 2-prong fitter (only once)
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setWeightedFinalPCA(false);
    fitter.setBz(magneticField);
    fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrNONE);

    if (doprocessFindDmesons) {
      histos.add("h2dGenD", "h2dGenD", kTH2F, {axisPt, axisEta});
      histos.add("h2dGenDbar", "h2dGenDbar", kTH2F, {axisPt, axisEta});
      histos.add("h3dRecD", "h2dRecD", kTH3F, {axisPt, axisEta, axisDMass});
      histos.add("h3dRecDbar", "h2dRecDbar", kTH3F, {axisPt, axisEta, axisDMass});

      histos.add("hMassD", "hMassD", kTH1F, {axisDMass});
      histos.add("hMassDbar", "hMassDbar", kTH1F, {axisDMass});

      if (doDCAplots) {
        histos.add("h2dDCAxyVsPtPiPlusFromD", "h2dDCAxyVsPtPiPlusFromD", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtPiMinusFromD", "h2dDCAxyVsPtPiMinusFromD", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtKaPlusFromD", "h2dDCAxyVsPtKaPlusFromD", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtKaMinusFromD", "h2dDCAxyVsPtKaMinusFromD", kTH2F, {axisPt, axisDCA});
      }
    }
  }

  void processFindDmesons(aod::Collision const& collision, alice3tracks const& tracks, aod::McParticles const& mcParticles)
  {
    // no grouping for MC particles -> as intended
    for (auto const& mcParticle : trueD)
      histos.fill(HIST("h2dGenD"), mcParticle.pt(), mcParticle.eta());
    for (auto const& mcParticle : trueDbar)
      histos.fill(HIST("h2dGenDbar"), mcParticle.pt(), mcParticle.eta());

    // group with this collision
    auto tracksPiPlusFromDgrouped = tracksPiPlusFromD->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksKaMinusFromDgrouped = tracksKaMinusFromD->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksKaPlusFromDgrouped = tracksKaPlusFromD->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksPiMinusFromDgrouped = tracksPiMinusFromD->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // D mesons
    for (auto const& posTrackRow : tracksPiPlusFromDgrouped) {
      for (auto const& negTrackRow : tracksKaMinusFromDgrouped) {
        if (mcSameMotherCheck && !checkSameMother(posTrackRow, negTrackRow))
          continue;
        if (!buildDecayCandidate(posTrackRow, negTrackRow, o2::constants::physics::MassPionCharged, o2::constants::physics::MassKaonCharged))
          continue;
        histos.fill(HIST("hMassD"), dmeson.mass);
        histos.fill(HIST("h3dRecD"), dmeson.pt, dmeson.eta, dmeson.mass);
      }
    }
    // D mesons
    for (auto const& posTrackRow : tracksKaPlusFromDgrouped) {
      for (auto const& negTrackRow : tracksPiMinusFromDgrouped) {
        if (mcSameMotherCheck && !checkSameMother(posTrackRow, negTrackRow))
          continue;
        if (!buildDecayCandidate(posTrackRow, negTrackRow, o2::constants::physics::MassKaonCharged, o2::constants::physics::MassPionCharged))
          continue;
        histos.fill(HIST("hMassDbar"), dmeson.mass);
        histos.fill(HIST("h3dRecDbar"), dmeson.pt, dmeson.eta, dmeson.mass);
      }
    }
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
  PROCESS_SWITCH(alice3decayFinder, processFindDmesons, "find D mesons", true);
  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<alice3decayPreselector>(cfgc),
    adaptAnalysisTask<alice3decayFinder>(cfgc)};
}
