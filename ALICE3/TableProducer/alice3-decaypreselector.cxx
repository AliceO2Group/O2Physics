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
#include <vector>
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
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/A3DecayFinderTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// simple checkers
// #define biton(var, nbit) ((var) |= (static_cast<uint32_t>(1) << (nbit)))
#define bitoff(var, nbit) ((var) &= ~(static_cast<uint32_t>(1) << (nbit))) //((a) &= ~(1ULL<<(b)))
// #define bitcheck(var, nbit) ((var) & (static_cast<uint32_t>(1) << (nbit)))

using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov>;

// For MC association in pre-selection
using labeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using tofTracks = soa::Join<aod::Tracks, aod::UpgradeTofs>;
using richTracks = soa::Join<aod::Tracks, aod::UpgradeRichs, aod::UpgradeRichSignal>;

struct alice3decaypreselector {
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

  void init(InitContext&)
  {
    // future dev if needed
  }

  // go declarative: use partitions instead of "if", then just toggle bits to allow for mask selection later
  Partition<tofTracks> pInnerTOFPi = nabs(aod::upgrade_tof::nSigmaPionInnerTOF) > nSigmaTOF;
  Partition<tofTracks> pInnerTOFKa = nabs(aod::upgrade_tof::nSigmaKaonInnerTOF) > nSigmaTOF;
  Partition<tofTracks> pInnerTOFPr = nabs(aod::upgrade_tof::nSigmaProtonInnerTOF) > nSigmaTOF;
  Partition<tofTracks> pOuterTOFPi = nabs(aod::upgrade_tof::nSigmaPionOuterTOF) > nSigmaTOF;
  Partition<tofTracks> pOuterTOFKa = nabs(aod::upgrade_tof::nSigmaKaonOuterTOF) > nSigmaTOF;
  Partition<tofTracks> pOuterTOFPr = nabs(aod::upgrade_tof::nSigmaProtonOuterTOF) > nSigmaTOF;
  Partition<richTracks> pRICHPi = aod::upgrade_rich::hasSig && aod::upgrade_rich::hasSigPi && nabs(aod::upgrade_rich::nSigmaPionRich) > nSigmaRICH;
  Partition<richTracks> pRICHKa = aod::upgrade_rich::hasSig && aod::upgrade_rich::hasSigKa && nabs(aod::upgrade_rich::nSigmaKaonRich) > nSigmaRICH;
  Partition<richTracks> pRICHPr = aod::upgrade_rich::hasSig && aod::upgrade_rich::hasSigPr && nabs(aod::upgrade_rich::nSigmaProtonRich) > nSigmaRICH;

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
  void processFilterInnerTOF(tofTracks const&)
  {
    for (auto const& track : pInnerTOFPi)
      bitoff(selectionMap[track.globalIndex()], kInnerTOFPion);
    for (auto const& track : pInnerTOFKa)
      bitoff(selectionMap[track.globalIndex()], kInnerTOFKaon);
    for (auto const& track : pInnerTOFPr)
      bitoff(selectionMap[track.globalIndex()], kInnerTOFProton);
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFilterOuterTOF(tofTracks const&)
  {
    for (auto const& track : pOuterTOFPi)
      bitoff(selectionMap[track.globalIndex()], kOuterTOFPion);
    for (auto const& track : pOuterTOFKa)
      bitoff(selectionMap[track.globalIndex()], kOuterTOFKaon);
    for (auto const& track : pOuterTOFPr)
      bitoff(selectionMap[track.globalIndex()], kOuterTOFProton);
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFilterRICH(richTracks const&)
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
      // D mesons
      if (!checkPDG(track, 421, -321)) //+421 -> -321 +211
        bitoff(selectionMap[track.globalIndex()], kTrueKaMinusFromD);
      if (!checkPDG(track, -421, +321)) //-421 -> +321 -211
        bitoff(selectionMap[track.globalIndex()], kTrueKaPlusFromD);
      if (!checkPDG(track, 421, +211)) //+421 -> -321 +211
        bitoff(selectionMap[track.globalIndex()], kTruePiPlusFromD);
      if (!checkPDG(track, -421, -211)) //-421 -> +321 -211
        bitoff(selectionMap[track.globalIndex()], kTruePiMinusFromD);

      // Lambdac baryons
      if (!checkPDG(track, +4122, +2212)) //+4122 -> +2212 -321 +211
        bitoff(selectionMap[track.globalIndex()], kTruePrPlusFromLc);
      if (!checkPDG(track, +4122, -321)) //+4122 -> +2212 -321 +211
        bitoff(selectionMap[track.globalIndex()], kTrueKaMinusFromLc);
      if (!checkPDG(track, +4122, +211)) //+4122 -> +2212 -321 +211
        bitoff(selectionMap[track.globalIndex()], kTruePiPlusFromLc);
      if (!checkPDG(track, -4122, -2212)) //-4122 -> -2212 +321 -211
        bitoff(selectionMap[track.globalIndex()], kTruePrMinusFromLc);
      if (!checkPDG(track, -4122, +321)) //-4122 -> -2212 +321 -211
        bitoff(selectionMap[track.globalIndex()], kTrueKaPlusFromLc);
      if (!checkPDG(track, -4122, -211)) //-4122 -> -2212 +321 -211
        bitoff(selectionMap[track.globalIndex()], kTruePiMinusFromLc);

      // XiCC daughters
      if (!checkPDG(track, 4422, 211)) // 4422 -> 4232 211, pi from xicc
        bitoff(selectionMap[track.globalIndex()], kTruePiFromXiCC);
      if (!checkPDG(track, 4232, 3312)) // 4232 -> 3312 211 211, xi from xic
        bitoff(selectionMap[track.globalIndex()], kTrueXiFromXiC);
      if (!checkPDG(track, 4232, 211)) // 4232 -> 3312 211 211, pi from xic
        bitoff(selectionMap[track.globalIndex()], kTruePiFromXiC);
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
  PROCESS_SWITCH(alice3decaypreselector, processInitialize, "Initialize (MUST be on)", true);
  PROCESS_SWITCH(alice3decaypreselector, processFilterInnerTOF, "Switch to use inner TOF PID", false);
  PROCESS_SWITCH(alice3decaypreselector, processFilterOuterTOF, "Switch to use outer TOF PID", false);
  PROCESS_SWITCH(alice3decaypreselector, processFilterRICH, "Switch to use RICH", false);
  PROCESS_SWITCH(alice3decaypreselector, processFilterOnMonteCarloTruth, "Switch to use MC truth", false);
  PROCESS_SWITCH(alice3decaypreselector, processPublishDecision, "Fill decision mask table (MUST be on)", true);
  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<alice3decaypreselector>(cfgc)};
}
