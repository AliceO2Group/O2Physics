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

/// \file alice3MulticharmFinder.cxx
/// \brief produces table of xicc candidates
/// \author Jesper Karlsson Gumprecht <jesper.gumprecht@cern.ch>

//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//   Decay finder task for ALICE 3
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Uses specific ALICE 3 PID and performance for studying
//    HF decays. Work in progress: use at your own risk!

#include "ALICE3/DataModel/A3DecayFinderTables.h"
#include "ALICE3/DataModel/OTFCollision.h"
#include "ALICE3/DataModel/OTFMulticharm.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsVertexing/PVertexer.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// simple checkers
// #define biton(var, nbit) ((var) |= (static_cast<uint32_t>(1) << (nbit)))
// #define bitoff(var, nbit) ((var) &= ~(static_cast<uint32_t>(1) << (nbit))) //((a) &= ~(1ULL<<(b)))
#define BIT_CHECK(var, nbit) ((var) & (static_cast<uint32_t>(1) << (nbit)))
#define GET_HIST(type, name) std::get<std::shared_ptr<type>>(histPointers[name])

using Alice3Tracks = soa::Join<aod::Tracks, aod::TracksCov, aod::Alice3DecayMaps, aod::McTrackLabels, aod::TracksDCA, aod::TracksExtraA3, aod::UpgradeTofs, aod::UpgradeTofExpectedTimes, aod::UpgradeRichs, aod::UpgradeRichSignals>;
using Alice3Collision = soa::Join<aod::Collisions, aod::OTFLUTConfigId>;

struct Alice3MulticharmFinder {
  SliceCache cache;

  Produces<aod::MCharmIndices> multiCharmIdx;
  Produces<aod::MCharmCores> multiCharmCore;
  Produces<aod::MCharmPID> multiCharmPID;
  Produces<aod::MCharmExtra> multiCharmExtra;

  // Operation and minimisation criteria
  Configurable<bool> fillDerivedTable{"fillDerivedTable", false, "fill MCharm[] tables (careful: memory)"};
  Configurable<float> magneticField{"magneticField", 20.0f, "Magnetic field (in kilogauss)"};
  Configurable<bool> doDCAplots{"doDCAplots", true, "do daughter prong DCA plots for D mesons"};
  Configurable<bool> mcSameMotherCheck{"mcSameMotherCheck", true, "check if tracks come from the same MC mother"};
  Configurable<std::vector<float>> minNTracks{"minNTracks", {-1}, "Minimum number of tracks"};

  Configurable<float> xiMinConstDCAxy{"xiMinConstDCAxy", 0.0005f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinConstDCAz{"xiMinConstDCAz", 0.0005f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinPtDepDCAxy{"xiMinPtDepDCAxy", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinPtDepDCAz{"xiMinPtDepDCAz", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinDecayRadius{"xiMinDecayRadius", 0.5, "Minimum R2D for XiC decay (cm)"};
  Configurable<float> xiMassWindow{"xiMassWindow", 0.005, "Mass window around Xi peak (GeV/c)"};

  Configurable<float> picTofDiffInner{"picTofDiffInner", 99999, "|signal - expected| (ps)"};
  Configurable<float> picMinConstDCAxy{"picMinConstDCAxy", 0.0005f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> picMinConstDCAz{"picMinConstDCAz", 0.0005f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> picMinPtDepDCAxy{"picMinPtDepDCAxy", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> picMinPtDepDCAz{"picMinPtDepDCAz", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> picMinPt{"picMinPt", 0.15, "Minimum pT for XiC pions"};

  Configurable<float> xicMaxDauDCA{"xicMaxDauDCA", 0.005f, "DCA between XiC daughters (cm)"};
  Configurable<float> xicMaxDCAxy{"xicMaxDCAxy", 0.0005f, "maxDCA"};
  Configurable<float> xicMaxDCAz{"xicMaxDCAz", 0.0005f, "maxDCA"};
  Configurable<float> xicMinDecayRadius{"xicMinDecayRadius", -1, "Minimum R2D for XiC decay (cm)"};
  Configurable<float> xicMinDecayDistanceFromPV{"xicMinDecayDistanceFromPV", -1, "Minimum distance for XiC decay from PV (cm)"};
  Configurable<float> xicMinProperLength{"xicMinProperLength", 0.002, "Minimum proper length for XiC decay (cm)"};
  Configurable<float> xicMaxProperLength{"xicMaxProperLength", 0.1, "Minimum proper length for XiC decay (cm)"};
  Configurable<float> xicMassWindow{"xicMassWindow", 0.012, "Mass window around XiC peak (GeV/c)"};

  Configurable<float> piccTofDiffInner{"piccTofDiffInner", 99999, "|signal - expected| (ps)"};
  Configurable<float> piccMinConstDCAxy{"piccMinConstDCAxy", 0.0005f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piccMinConstDCAz{"piccMinConstDCAz", 0.0005f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piccMinPtDepDCAxy{"piccMinPtDepDCAxy", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piccMinPtDepDCAz{"piccMinPtDepDCAz", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piccMinPt{"piccMinPt", 0.3, "Minimum pT for XiCC pions"};

  Configurable<float> xiccMaxDauDCA{"xiccMaxDauDCA", 0.005f, "DCA between XiCC daughters (cm)"};
  Configurable<float> xiccMaxDCAxy{"xiccMaxDCAxy", 0.005f, "maxDCA"};
  Configurable<float> xiccMaxDCAz{"xiccMaxDCAz", 0.005f, "maxDCA"};
  Configurable<float> xiccMaxEta{"xiccMaxEta", 1.5, "Max eta"};
  Configurable<float> xiccMinDecayRadius{"xiccMinDecayRadius", -1, "Minimum R2D for XiCC decay (cm)"};
  Configurable<float> xiccMinProperLength{"xiccMinProperLength", -1, "Minimum proper length for XiCC decay (cm)"};
  Configurable<float> xiccMaxProperLength{"xiccMaxProperLength", 999, "Minimum proper length for XiCC decay (cm)"};
  Configurable<float> xiccMassWindow{"xiccMassWindow", 0.25, "Mass window around XiCC peak (GeV/c). Make sure that bkg region is included in this window"};

  ConfigurableAxis axisEta{"axisEta", {80, -4.0f, +4.0f}, "#eta"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisDCA2D{"axisDCA2D", {400, -200, 200}, "DCA2d (#mum)"};
  ConfigurableAxis axisDCA{"axisDCA", {400, 0, 400}, "DCA (#mum)"};
  ConfigurableAxis axisRadius{"axisRadius", {10000, 0, 10000}, "Decay radius (#mum)"};
  ConfigurableAxis axisRadius2D{"axisRadius2D", {1000, 0, 100000}, "Decay radius (#mum)"};
  ConfigurableAxis axisRadius2DXi{"axisRadius2DXi", {1000, 0, 20}, "Decay radius (cm)"};
  ConfigurableAxis axisDecayLength{"axisDecayLength", {2000, 0, 2000}, "Decay lenght (#mum)"};
  ConfigurableAxis axisTOFTrack{"axisTOFTrack", {1000, 0, 5000}, "TOF track time"};

  ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.221f, 1.421f}, "Xi Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisXiCMass{"axisXiCMass", {200, 2.368f, 2.568f}, "XiC Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisXiCCMass{"axisXiCCMass", {200, 3.521f, 3.721f}, "XiCC Inv Mass (GeV/c^{2})"};

  ConfigurableAxis axisDCAXiCDaughters{"axisDCAXiCDaughters", {200, 0, 100}, "DCA (mum)"};
  ConfigurableAxis axisDCAXiCCDaughters{"axisDCAXiCCDaughters", {200, 0, 100}, "DCA (mum)"};

  ConfigurableAxis axisNConsidered{"axisNConsidered", {200, -0.5f, 199.5f}, "Number of considered track combinations"};

  o2::vertexing::DCAFitterN<2> fitter;
  o2::vertexing::DCAFitterN<3> fitter3;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::map<std::string, HistPtr> histPointers;
  std::vector<int> savedConfigs;

  // Constants
  static constexpr std::array<int, 6> MomentumIndices = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component

  // filter expressions for pions
  static constexpr uint32_t TrackSelectionPic = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiFromXiC;
  static constexpr uint32_t TrackSelectionPicc = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiFromXiCC;

  // partitions
  Partition<aod::McParticles> trueXi = aod::mcparticle::pdgCode == o2::constants::physics::MassXiMinus;
  Partition<aod::McParticles> trueXiC = aod::mcparticle::pdgCode == o2::constants::physics::MassXiCPlus;
  Partition<aod::McParticles> trueXiCC = aod::mcparticle::pdgCode == o2::constants::physics::MassXiCCPlusPlus;

  Partition<Alice3Tracks> picTracks =
    ((aod::a3DecayMap::decayMap & TrackSelectionPic) == TrackSelectionPic) && aod::track::signed1Pt > 0.0f && 1.0f / nabs(aod::track::signed1Pt) > picMinPt&& nabs(aod::track::dcaXY) > picMinConstDCAxy + picMinPtDepDCAxy* nabs(aod::track::signed1Pt) && nabs(aod::track::dcaZ) > picMinConstDCAz + picMinPtDepDCAz* nabs(aod::track::signed1Pt);

  Partition<Alice3Tracks> piccTracks =
    ((aod::a3DecayMap::decayMap & TrackSelectionPicc) == TrackSelectionPicc) && aod::track::signed1Pt > 0.0f && 1.0f / nabs(aod::track::signed1Pt) > piccMinPt&& nabs(aod::track::dcaXY) > piccMinConstDCAxy + piccMinPtDepDCAxy* nabs(aod::track::signed1Pt) && nabs(aod::track::dcaZ) > piccMinConstDCAz + piccMinPtDepDCAz* nabs(aod::track::signed1Pt);

  // Helper struct to pass candidate information
  struct {
    // decay properties
    float dca;
    float mass;
    float pt;
    float eta;
    std::array<float, 3> xyz;
    std::array<float, 3> prong0mom;
    std::array<float, 3> prong1mom;
    std::array<float, 3> prong2mom;
    std::array<float, o2::track::kLabCovMatSize> parentTrackCovMatrix;
  } thisXiCcandidate;

  struct {
    float dca;
    float mass;
    float pt;
    float eta;
    std::array<float, 3> xyz;
    std::array<float, 3> prong0mom;
    std::array<float, 3> prong1mom;
    std::array<float, o2::track::kLabCovMatSize> parentTrackCovMatrix;

    // charm daughters
    int nSiliconHitsPiCC;
    int nTPCHitsPiCC;
  } thisXiCCcandidate;

  template <typename TTrackType>
  bool buildDecayCandidateTwoBody(TTrackType const& t0, TTrackType const& t1, float mass0, float mass1)
  {
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter.process(t0, t1);
    } catch (...) {
      return false;
    }
    if (nCand == 0) {
      return false;
    }
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}

    o2::track::TrackParCov t0new = fitter.getTrack(0);
    o2::track::TrackParCov t1new = fitter.getTrack(1);
    t0new.getPxPyPzGlo(thisXiCCcandidate.prong0mom);
    t1new.getPxPyPzGlo(thisXiCCcandidate.prong1mom);

    // get decay vertex coordinates
    const auto& vtx = fitter.getPCACandidate();
    for (size_t i = 0; i < thisXiCCcandidate.xyz.size(); i++) {
      thisXiCCcandidate.xyz[i] = vtx[i];
    }

    // compute cov mat
    for (int ii = 0; ii < o2::track::kLabCovMatSize; ii++)
      thisXiCCcandidate.parentTrackCovMatrix[ii] = 0.0f;

    std::array<float, o2::track::kLabCovMatSize> covA = {0};
    std::array<float, o2::track::kLabCovMatSize> covB = {0};
    fitter.getTrack(0).getCovXYZPxPyPzGlo(covA);
    fitter.getTrack(1).getCovXYZPxPyPzGlo(covB);

    for (size_t i = 0; i < MomentumIndices.size(); i++) {
      int j = MomentumIndices[i];
      thisXiCCcandidate.parentTrackCovMatrix[j] = covA[j] + covB[j];
    }

    auto covVtx = fitter.calcPCACovMatrix();
    thisXiCCcandidate.parentTrackCovMatrix[0] = covVtx(0, 0);
    thisXiCCcandidate.parentTrackCovMatrix[1] = covVtx(1, 0);
    thisXiCCcandidate.parentTrackCovMatrix[2] = covVtx(1, 1);
    thisXiCCcandidate.parentTrackCovMatrix[3] = covVtx(2, 0);
    thisXiCCcandidate.parentTrackCovMatrix[4] = covVtx(2, 1);
    thisXiCCcandidate.parentTrackCovMatrix[5] = covVtx(2, 2);

    // set relevant values
    thisXiCCcandidate.dca = std::sqrt(fitter.getChi2AtPCACandidate());
    if (thisXiCCcandidate.dca > xiccMaxDauDCA) {
      return false;
    }

    thisXiCCcandidate.mass = RecoDecay::m(std::array{std::array{thisXiCCcandidate.prong0mom[0], thisXiCCcandidate.prong0mom[1], thisXiCCcandidate.prong0mom[2]}, std::array{thisXiCCcandidate.prong1mom[0], thisXiCCcandidate.prong1mom[1], thisXiCCcandidate.prong1mom[2]}}, std::array{mass0, mass1});

    if (std::fabs(thisXiCCcandidate.mass - o2::constants::physics::MassXiCCPlusPlus) > xiccMassWindow) {
      return false;
    }

    thisXiCCcandidate.pt = std::hypot(thisXiCCcandidate.prong0mom[0] + thisXiCCcandidate.prong1mom[0], thisXiCCcandidate.prong0mom[1] + thisXiCCcandidate.prong1mom[1]);
    thisXiCCcandidate.eta = RecoDecay::eta(std::array{thisXiCCcandidate.prong0mom[0] + thisXiCCcandidate.prong1mom[0], thisXiCCcandidate.prong0mom[1] + thisXiCCcandidate.prong1mom[1], thisXiCCcandidate.prong0mom[2] + thisXiCCcandidate.prong1mom[2]});
    return true;
  }

  template <typename TTrackType1, typename TTrackType2, typename TTrackType3>
  bool buildDecayCandidateThreeBody(TTrackType1 const& prong0, TTrackType2 const& prong1, TTrackType3 const& prong2, float p0mass, float p1mass, float p2mass)
  {
    o2::track::TrackParCov t0 = getTrackParCov(prong0);
    o2::track::TrackParCov t1 = getTrackParCov(prong1);
    o2::track::TrackParCov t2 = getTrackParCov(prong2);

    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter3.process(t0, t1, t2);
    } catch (...) {
      return false;
    }
    if (nCand == 0) {
      return false;
    }
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}

    t0 = fitter3.getTrack(0);
    t1 = fitter3.getTrack(1);
    t2 = fitter3.getTrack(2);
    t0.getPxPyPzGlo(thisXiCcandidate.prong0mom);
    t1.getPxPyPzGlo(thisXiCcandidate.prong1mom);
    t2.getPxPyPzGlo(thisXiCcandidate.prong2mom);

    // get decay vertex coordinates
    const auto& vtx = fitter3.getPCACandidate();
    for (size_t i = 0; i < thisXiCcandidate.xyz.size(); i++) {
      thisXiCcandidate.xyz[i] = vtx[i];
    }

    // compute cov mat
    for (int ii = 0; ii < o2::track::kLabCovMatSize; ii++)
      thisXiCcandidate.parentTrackCovMatrix[ii] = 0.0f;

    std::array<float, o2::track::kLabCovMatSize> covA = {0};
    std::array<float, o2::track::kLabCovMatSize> covB = {0};
    std::array<float, o2::track::kLabCovMatSize> covC = {0};
    fitter3.getTrack(0).getCovXYZPxPyPzGlo(covA);
    fitter3.getTrack(1).getCovXYZPxPyPzGlo(covB);
    fitter3.getTrack(2).getCovXYZPxPyPzGlo(covC);

    for (size_t i = 0; i < MomentumIndices.size(); i++) {
      int j = MomentumIndices[i];
      thisXiCcandidate.parentTrackCovMatrix[j] = covA[j] + covB[j] + covC[j];
    }

    auto covVtx = fitter3.calcPCACovMatrix();
    thisXiCcandidate.parentTrackCovMatrix[0] = covVtx(0, 0);
    thisXiCcandidate.parentTrackCovMatrix[1] = covVtx(1, 0);
    thisXiCcandidate.parentTrackCovMatrix[2] = covVtx(1, 1);
    thisXiCcandidate.parentTrackCovMatrix[3] = covVtx(2, 0);
    thisXiCcandidate.parentTrackCovMatrix[4] = covVtx(2, 1);
    thisXiCcandidate.parentTrackCovMatrix[5] = covVtx(2, 2);

    // set relevant values
    thisXiCcandidate.dca = std::sqrt(fitter3.getChi2AtPCACandidate());
    if (thisXiCcandidate.dca > xicMaxDauDCA) {
      return false;
    }
    thisXiCcandidate.mass = RecoDecay::m(std::array{std::array{thisXiCcandidate.prong0mom[0], thisXiCcandidate.prong0mom[1], thisXiCcandidate.prong0mom[2]}, std::array{thisXiCcandidate.prong1mom[0], thisXiCcandidate.prong1mom[1], thisXiCcandidate.prong1mom[2]}, std::array{thisXiCcandidate.prong2mom[0], thisXiCcandidate.prong2mom[1], thisXiCcandidate.prong2mom[2]}}, std::array{p0mass, p1mass, p2mass});
    thisXiCcandidate.pt = std::hypot(thisXiCcandidate.prong0mom[0] + thisXiCcandidate.prong1mom[0] + thisXiCcandidate.prong2mom[0], thisXiCcandidate.prong0mom[1] + thisXiCcandidate.prong1mom[1] + thisXiCcandidate.prong2mom[1]);
    thisXiCcandidate.eta = RecoDecay::eta(std::array{thisXiCcandidate.prong0mom[0] + thisXiCcandidate.prong1mom[0] + thisXiCcandidate.prong2mom[0], thisXiCcandidate.prong0mom[1] + thisXiCcandidate.prong1mom[1] + thisXiCcandidate.prong2mom[1], thisXiCcandidate.prong0mom[2] + thisXiCcandidate.prong1mom[2] + thisXiCcandidate.prong2mom[2]});
    return true;
  }

  template <typename TTrackType>
  int getPdgCodeForTrack(const TTrackType& track)
  {
    auto mcParticle = track.template mcParticle_as<aod::McParticles>();
    return mcParticle.pdgCode();
  }

  /// function to check if tracks have the same mother in MC
  template <typename TTrackType1, typename TTrackType2>
  bool checkSameMother(TTrackType1 const& track1, TTrackType2 const& track2)
  {
    bool returnValue = false;
    // Association check
    // There might be smarter ways of doing this in the future
    if (track1.has_mcParticle() && track2.has_mcParticle()) {
      auto mcParticle1 = track1.template mcParticle_as<aod::McParticles>();
      auto mcParticle2 = track2.template mcParticle_as<aod::McParticles>();
      if (mcParticle1.has_mothers() && mcParticle2.has_mothers()) {
        for (const auto& mcParticleMother1 : mcParticle1.template mothers_as<aod::McParticles>()) {
          for (const auto& mcParticleMother2 : mcParticle2.template mothers_as<aod::McParticles>()) {
            if (mcParticleMother1.globalIndex() == mcParticleMother2.globalIndex()) {
              returnValue = true;
            }
          }
        }
      }
    } // end association check
    return returnValue;
  }

  // Association check for the XiCC pion
  template <typename TTrackType1, typename TTrackType2>
  bool checkSameMotherExtra(TTrackType1 const& track1, TTrackType2 const& track2)
  {
    bool returnValue = false;
    // This might perhaps be a bit excessive
    // Could be joined with `checkSameMother` but leaving as is for now
    if (track1.has_mcParticle() && track2.has_mcParticle()) {
      auto mcParticle1 = track1.template mcParticle_as<aod::McParticles>();
      auto mcParticle2 = track2.template mcParticle_as<aod::McParticles>();
      if (mcParticle1.has_mothers() && mcParticle2.has_mothers()) {
        for (const auto& mcParticleMother1 : mcParticle1.template mothers_as<aod::McParticles>()) {
          if (mcParticleMother1.has_mothers()) {
            for (const auto& mcParticleGrandMother1 : mcParticleMother1.template mothers_as<aod::McParticles>()) {
              for (const auto& mcParticleMother2 : mcParticle2.template mothers_as<aod::McParticles>()) {
                if (mcParticleGrandMother1.globalIndex() == mcParticleMother2.globalIndex()) {
                  returnValue = true;
                }
              }
            }
          }
        }
      }
    } // end association check
    return returnValue;
  }

  void init(InitContext&)
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

    fitter3.setPropagateToPCA(true);
    fitter3.setMaxR(200.);
    fitter3.setMinParamChange(1e-3);
    fitter3.setMinRelChi2Change(0.9);
    fitter3.setMaxDZIni(1e9);
    fitter3.setMaxChi2(1e9);
    fitter3.setUseAbsDCA(true);
    fitter3.setWeightedFinalPCA(false);
    fitter3.setBz(magneticField);
    fitter3.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrNONE);

    // This histogram bookkeeps the attempts at DCA minimization and their eventual
    // failure rates.
    // --- 0: attempt XiC, 1: success XiC
    // --- 2: attempt XiCC, 3: success XiCC
    histos.add("hCharmBuilding", "hCharmBuilding", kTH1D, {{10, -0.5, 9.5f}});

    histos.add("h2dGenXi", "h2dGenXi", kTH2D, {axisPt, axisEta});
    histos.add("h2dGenXiC", "h2dGenXiC", kTH2D, {axisPt, axisEta});
    histos.add("h2dGenXiCC", "h2dGenXiCC", kTH2D, {axisPt, axisEta});

    histos.add("hMassXi", "hMassXi", kTH1D, {axisXiMass});
    histos.add("hMassXiC", "hMassXiC", kTH1D, {axisXiCMass});
    histos.add("hMassXiCC", "hMassXiCC", kTH1D, {axisXiCCMass});

    histos.add("hEtaXiCC", "hEtaXiCC", kTH1D, {axisEta});
    histos.add("hPtXiCC", "hPtXiCC", kTH1D, {axisPt});
    histos.add("h3dMassXiCC", "h3dMassXiCC", kTH3D, {axisPt, axisEta, axisXiCCMass});

    histos.add("hDCAXiCDaughters", "hDCAXiCDaughters", kTH1D, {axisDCAXiCDaughters});
    histos.add("hDCAXiCCDaughters", "hDCAXiCCDaughters", kTH1D, {axisDCAXiCCDaughters});
    histos.add("hDCAxyXi", "hDCAxyXi", kTH1D, {axisDCA});
    histos.add("hDCAzXi", "hDCAzXi", kTH1D, {axisDCA});

    histos.add("hDCAxyXiC", "hDCAxyXiC", kTH1D, {axisDCA});
    histos.add("hDCAzXiC", "hDCAzXiC", kTH1D, {axisDCA});

    histos.add("hDCAxyXiCC", "hDCAxyXiCC", kTH1D, {axisDCA});
    histos.add("hDCAzXiCC", "hDCAzXiCC", kTH1D, {axisDCA});

    histos.add("hPi1cPt", "hPi1cPt", kTH1D, {axisPt});
    histos.add("hPi2cPt", "hPi2cPt", kTH1D, {axisPt});
    histos.add("hPiccPt", "hPiccPt", kTH1D, {axisPt});

    histos.add("hPi1cDCAxy", "hPi1cDCAxy", kTH1D, {axisDCA});
    histos.add("hPi1cDCAz", "hPi1cDCAz", kTH1D, {axisDCA});
    histos.add("hPi2cDCAxy", "hPi2cDCAxy", kTH1D, {axisDCA});
    histos.add("hPi2cDCAz", "hPi2cDCAz", kTH1D, {axisDCA});
    histos.add("hPiccDCAxy", "hPiccDCAxy", kTH1D, {axisDCA});
    histos.add("hPiccDCAz", "hPiccDCAz", kTH1D, {axisDCA});

    histos.add("hMinXiDecayRadius", "hMinXiDecayRadius", kTH1D, {axisRadius2DXi});
    histos.add("hMinXiCDecayRadius", "hMinXiCDecayRadius", kTH1D, {axisRadius});
    histos.add("hMinXiCCDecayRadius", "hMinXiCCDecayRadius", kTH1D, {axisRadius});

    histos.add("hMinxicDecayDistanceFromPV", "hMinxicDecayDistanceFromPV", kTH1D, {axisDecayLength});
    histos.add("hProperLengthXiC", "hProperLengthXiC", kTH1D, {axisDecayLength});
    histos.add("hProperLengthXiCC", "hProperLengthXiCC", kTH1D, {axisDecayLength});

    histos.add("hInnerTOFTrackTimeRecoPi1c", "hInnerTOFTrackTimeRecoPi1c", kTH1D, {axisTOFTrack});
    histos.add("hInnerTOFTrackTimeRecoPi2c", "hInnerTOFTrackTimeRecoPi2c", kTH1D, {axisTOFTrack});
    histos.add("hInnerTOFTrackTimeRecoPicc", "hInnerTOFTrackTimeRecoPicc", kTH1D, {axisTOFTrack});

    histos.add("hXiRadiusVsXicRadius", "hXiRadiusVsXicRadius", kTH2D, {axisRadius2D, axisRadius2D});
    histos.add("hXicRadiusVsXiccRadius", "hXicRadiusVsXiccRadius", kTH2D, {axisRadius2D, axisRadius2D});

    // These histograms bookkeep the exact number of combinations attempted
    // CombinationsXiC: triplets Xi-pi-pi considered per Xi
    // CombinationsXiCC: doublets XiC-pi considered per XiC
    histos.add("hCombinationsXiC", "hCombinationsXiC", kTH1D, {axisNConsidered});
    histos.add("hCombinationsXiCC", "hCombinationsXiCC", kTH1D, {axisNConsidered});
    histos.add("hNCollisions", "hNCollisions", kTH1D, {{2, 0.5, 2.5}});
    histos.add("hNTracks", "hNTracks", kTH1D, {{20000, 0, 20000}});

    if (doDCAplots) {
      histos.add("h2dDCAxyVsPtXiFromXiC", "h2dDCAxyVsPtXiFromXiC", kTH2D, {axisPt, axisDCA2D});
      histos.add("h2dDCAxyVsPtPiFromXiC", "h2dDCAxyVsPtPiFromXiC", kTH2D, {axisPt, axisDCA2D});
      histos.add("h2dDCAxyVsPtPiFromXiCC", "h2dDCAxyVsPtPiFromXiCC", kTH2D, {axisPt, axisDCA2D});
    }
  }

  void initDetectorConfiguration(const int icfg)
  {
    if (std::find(savedConfigs.begin(), savedConfigs.end(), icfg) != savedConfigs.end()) {
      return;
    }

    savedConfigs.push_back(icfg);

    // do more plots
    const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";
    histPointers.insert({histPath + "hMassXiCC", histos.add((histPath + "hMassXiCC").c_str(), "hMassXiCC", {kTH1D, {{axisXiCCMass}}})});
    histPointers.insert({histPath + "hNCollisions", histos.add((histPath + "hNCollisions").c_str(), "hNCollisions", {kTH1D, {{2, 0.5, 2.5}}})});
    histPointers.insert({histPath + "hNTracks", histos.add((histPath + "hNTracks").c_str(), "hNTracks", {kTH1D, {{20000, 0, 20000}}})});
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processGenerated(aod::McParticles const&)
  {
    for (auto const& mcParticle : trueXi) {
      histos.fill(HIST("h2dGenXi"), mcParticle.pt(), mcParticle.eta());
    }
    for (auto const& mcParticle : trueXiC) {
      histos.fill(HIST("h2dGenXiC"), mcParticle.pt(), mcParticle.eta());
    }
    for (auto const& mcParticle : trueXiCC) {
      histos.fill(HIST("h2dGenXiCC"), mcParticle.pt(), mcParticle.eta());
    }
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFindXiCC(Alice3Collision::iterator const& collision, Alice3Tracks const& tracks, aod::McParticles const&, aod::UpgradeCascades const& cascades)
  {
    const std::string histPath = "Configuration_" + std::to_string(collision.lutConfigId()) + "/";
    initDetectorConfiguration(collision.lutConfigId());

    GET_HIST(TH1, histPath + "hNCollisions")->Fill(1);
    GET_HIST(TH1, histPath + "hNTracks")->Fill(tracks.size());
    std::cout << tracks.size() << std::endl;
    if (tracks.size() < minNTracks.value[collision.lutConfigId()]) {
      return;
    }
    GET_HIST(TH1, histPath + "hNCollisions")->Fill(2);

    // group with this collision
    // n.b. cascades do not need to be grouped, being used directly in iterator-grouping
    auto picTracksGrouped = picTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto piccTracksGrouped = piccTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);


    for (auto const& track : tracks) {
      if (BIT_CHECK(track.decayMap(), kTruePiFromXiC)) {
        histos.fill(HIST("h2dDCAxyVsPtPiFromXiC"), track.pt(), track.dcaXY() * 1e+4);
      }
      if (BIT_CHECK(track.decayMap(), kTruePiFromXiCC)) {
        histos.fill(HIST("h2dDCAxyVsPtPiFromXiCC"), track.pt(), track.dcaXY() * 1e+4);
      }
    }

    for (auto const& xiCand : cascades) {
      auto xi = xiCand.cascadeTrack_as<Alice3Tracks>(); // de-reference cascade track
      histos.fill(HIST("hMassXi"), xiCand.mXi());
      histos.fill(HIST("h2dDCAxyVsPtXiFromXiC"), xi.pt(), xi.dcaXY() * 1e+4);
      if (std::fabs(xiCand.mXi() - o2::constants::physics::MassXiMinus) > xiMassWindow) {
        continue; // out of mass region
      }

      uint32_t nCombinationsC = 0;
      auto bach = xiCand.bachTrack_as<Alice3Tracks>(); // de-reference bach track
      auto neg = xiCand.negTrack_as<Alice3Tracks>();   // de-reference neg track
      auto pos = xiCand.posTrack_as<Alice3Tracks>();   // de-reference pos track

      if (!BIT_CHECK(xi.decayMap(), kTrueXiFromXiC)) {
        continue;
      }

      if (std::fabs(xi.dcaXY()) < xiMinConstDCAxy || std::fabs(xi.dcaZ()) < xiMinConstDCAz) {
        continue; // likely a primary xi
      }

      histos.fill(HIST("hDCAxyXi"), std::fabs(xi.dcaXY() * 1e+4));
      histos.fill(HIST("hDCAzXi"), std::fabs(xi.dcaZ() * 1e+4));

      if (xiCand.cascRadius() < xiMinDecayRadius) {
        continue;
      }

      histos.fill(HIST("hMinXiDecayRadius"), xiCand.cascRadius());
      for (auto const& pi1c : picTracksGrouped) {
        if (mcSameMotherCheck && !checkSameMother(xi, pi1c)) {
          continue;
        }

        if (xiCand.posTrackId() == pi1c.globalIndex() || xiCand.negTrackId() == pi1c.globalIndex() || xiCand.bachTrackId() == pi1c.globalIndex()) {
          continue; // avoid using any track that was already used
        }

        if (pi1c.pt() < picMinPt) {
          continue; // too low momentum
        }

        histos.fill(HIST("hPi1cPt"), pi1c.pt());
        float pi1cTOFDiffInner = std::fabs(pi1c.innerTOFTrackTimeReco() - pi1c.innerTOFExpectedTimePi());
        float pi1cTOFDiffOuter = std::fabs(pi1c.outerTOFTrackTimeReco() - pi1c.outerTOFExpectedTimePi());
        if (pi1cTOFDiffInner > picTofDiffInner) {
          continue; // did not arrive at expected time
        }

        histos.fill(HIST("hInnerTOFTrackTimeRecoPi1c"), pi1cTOFDiffInner);
        // second pion from XiC decay for starts here
        for (auto const& pi2c : picTracksGrouped) {
          if (mcSameMotherCheck && !checkSameMother(xi, pi2c)) {
            continue; // keep only if same mother
          }

          if (pi1c.globalIndex() >= pi2c.globalIndex()) {
            continue; // avoid same-mother, avoid double-counting
          }

          if (xiCand.posTrackId() == pi2c.globalIndex() || xiCand.negTrackId() == pi2c.globalIndex() || xiCand.bachTrackId() == pi2c.globalIndex()) {
            continue; // avoid using any track that was already used
          }

          if (pi2c.pt() < picMinPt) {
            continue; // too low momentum
          }

          histos.fill(HIST("hPi2cPt"), pi2c.pt());
          float pi2cTOFDiffInner = std::fabs(pi2c.innerTOFTrackTimeReco() - pi2c.innerTOFExpectedTimePi());
          float pi2cTOFDiffOuter = std::fabs(pi2c.outerTOFTrackTimeReco() - pi2c.outerTOFExpectedTimePi());
          if (pi2cTOFDiffInner > picTofDiffInner) {
            continue; // did not arrive at expected time
          }

          histos.fill(HIST("hInnerTOFTrackTimeRecoPi2c"), pi2cTOFDiffInner);
          // if I am here, it means this is a triplet to be considered for XiC vertexing.
          // will now attempt to build a three-body decay candidate with these three track rows.

          nCombinationsC++;
          histos.fill(HIST("hCharmBuilding"), 0.0f);
          if (!buildDecayCandidateThreeBody(xi, pi1c, pi2c, o2::constants::physics::MassXiMinus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged)) {
            continue; // failed at building candidate
          }

          histos.fill(HIST("hDCAXiCDaughters"), thisXiCcandidate.dca * 1e+4);
          if (std::fabs(thisXiCcandidate.mass - o2::constants::physics::MassXiCPlus) > xicMassWindow) {
            continue; // out of mass region
          }

          histos.fill(HIST("hCharmBuilding"), 1.0f);

          const std::array<float, 3> momentumC = {
            thisXiCcandidate.prong0mom[0] + thisXiCcandidate.prong1mom[0] + thisXiCcandidate.prong2mom[0],
            thisXiCcandidate.prong0mom[1] + thisXiCcandidate.prong1mom[1] + thisXiCcandidate.prong2mom[1],
            thisXiCcandidate.prong0mom[2] + thisXiCcandidate.prong1mom[2] + thisXiCcandidate.prong2mom[2]};

          o2::track::TrackParCov xicTrack(thisXiCcandidate.xyz, momentumC, thisXiCcandidate.parentTrackCovMatrix, +1);
          float xicDecayRadius2D = std::hypot(thisXiCcandidate.xyz[0], thisXiCcandidate.xyz[1]);
          if (xicDecayRadius2D < xiccMinDecayRadius)
            continue; // do not take if radius too small, likely a primary combination

          histos.fill(HIST("hMinXiCDecayRadius"), xicDecayRadius2D * 1e+4);
          if (xicDecayRadius2D > xiCand.cascRadius()) {
            continue;
          }

          histos.fill(HIST("hXiRadiusVsXicRadius"), xiCand.cascRadius() * 1e+4, xicDecayRadius2D * 1e+4);
          o2::dataformats::DCA dcaInfo;
          float xicdcaXY = 1e+10, xicdcaZ = 1e+10;
          o2::track::TrackParCov xicTrackCopy(xicTrack); // paranoia
          o2::vertexing::PVertex primaryVertex;
          primaryVertex.setXYZ(collision.posX(), collision.posY(), collision.posZ());

          if (xicTrackCopy.propagateToDCA(primaryVertex, magneticField, &dcaInfo)) {
            xicdcaXY = dcaInfo.getY();
            xicdcaZ = dcaInfo.getZ();
          }

          if (std::fabs(xicdcaXY) < xicMaxDCAxy || std::fabs(xicdcaZ) < xicMaxDCAz) {
            continue; // likely a primary xic
          }

          histos.fill(HIST("hDCAxyXiC"), std::fabs(xicdcaXY * 1e+4));
          histos.fill(HIST("hDCAzXiC"), std::fabs(xicdcaZ * 1e+4));
          histos.fill(HIST("hMassXiC"), thisXiCcandidate.mass);

          // attempt XiCC finding
          uint32_t nCombinationsCC = 0;
          for (auto const& picc : piccTracksGrouped) {
            if (mcSameMotherCheck && !checkSameMotherExtra(xi, picc)) {
              continue;
            }

            if (xiCand.posTrackId() == picc.globalIndex() || xiCand.negTrackId() == picc.globalIndex() || xiCand.bachTrackId() == picc.globalIndex()) {
              continue; // avoid using any track that was already used
            }

            if (picc.pt() < piccMinPt) {
              continue; // too low momentum
            }

            histos.fill(HIST("hPiccPt"), picc.pt());
            float piccTOFDiffInner = std::fabs(picc.innerTOFTrackTimeReco() - picc.innerTOFExpectedTimePi());
            float piccTOFDiffOuter = std::fabs(picc.outerTOFTrackTimeReco() - picc.outerTOFExpectedTimePi());
            if (piccTOFDiffInner > piccTofDiffInner) {
              continue; // did not arrive at expected time
            }

            histos.fill(HIST("hInnerTOFTrackTimeRecoPicc"), piccTOFDiffInner);
            o2::track::TrackParCov piccTrack = getTrackParCov(picc);
            nCombinationsCC++;
            histos.fill(HIST("hCharmBuilding"), 2.0f);
            if (!buildDecayCandidateTwoBody(xicTrack, piccTrack, o2::constants::physics::MassXiCPlus, o2::constants::physics::MassPionCharged)) {
              continue; // failed at building candidate
            }

            histos.fill(HIST("hDCAXiCCDaughters"), thisXiCCcandidate.dca * 1e+4);

            const std::array<float, 3> momentumCC = {
              thisXiCCcandidate.prong0mom[0] + thisXiCCcandidate.prong1mom[0],
              thisXiCCcandidate.prong0mom[1] + thisXiCCcandidate.prong1mom[1],
              thisXiCCcandidate.prong0mom[2] + thisXiCCcandidate.prong1mom[2]};

            o2::track::TrackParCov xiccTrack(thisXiCCcandidate.xyz, momentumCC, thisXiCCcandidate.parentTrackCovMatrix, +2);
            float xiccDecayRadius2D = std::hypot(thisXiCCcandidate.xyz[0], thisXiCCcandidate.xyz[1]);
            if (xiccDecayRadius2D < xiccMinDecayRadius) {
              continue; // do not take if radius too small, likely a primary combination
            }

            histos.fill(HIST("hMinXiCCDecayRadius"), xiccDecayRadius2D * 1e+4);
            float totalMomentumC = std::hypot(momentumC[0], momentumC[1], momentumC[2]);
            float decayLengthXiC = std::hypot(
              thisXiCcandidate.xyz[0] - thisXiCCcandidate.xyz[0],
              thisXiCcandidate.xyz[1] - thisXiCCcandidate.xyz[1],
              thisXiCcandidate.xyz[2] - thisXiCCcandidate.xyz[2]);
            float xicProperLength = decayLengthXiC * thisXiCcandidate.mass / totalMomentumC;

            if (xicProperLength < xicMinProperLength || xicProperLength > xicMaxProperLength) {
              continue; // likely background
            }

            histos.fill(HIST("hProperLengthXiC"), xicProperLength * 1e+4);

            float xicDistanceFromPV = std::hypot(
              thisXiCcandidate.xyz[0] - collision.posX(),
              thisXiCcandidate.xyz[1] - collision.posY(),
              thisXiCcandidate.xyz[2] - collision.posZ());
            float xicDecayDistanceFromPV = xicDistanceFromPV * thisXiCcandidate.mass / totalMomentumC;
            if (xicDecayDistanceFromPV < xicMinDecayDistanceFromPV) {
              continue; // too close to PV
            }

            histos.fill(HIST("hMinxicDecayDistanceFromPV"), xicDecayDistanceFromPV * 1e+4);
            float totalMomentumCC = std::hypot(momentumCC[0], momentumCC[1], momentumCC[2]);
            float decayLengthXiCC = std::hypot(
              thisXiCCcandidate.xyz[0] - collision.posX(),
              thisXiCCcandidate.xyz[1] - collision.posY(),
              thisXiCCcandidate.xyz[2] - collision.posZ());
            float xiccProperLength = decayLengthXiCC * thisXiCCcandidate.mass / totalMomentumCC;
            if (xiccProperLength < xiccMinProperLength || xiccProperLength > xicMaxProperLength) {
              continue; // likely background
            }

            histos.fill(HIST("hProperLengthXiCC"), xiccProperLength * 1e+4);
            if (xiccDecayRadius2D > xicDecayRadius2D) {
              continue; // XiCC should decay before XiC
            }

            histos.fill(HIST("hXicRadiusVsXiccRadius"), xicDecayRadius2D * 1e+4, xiccDecayRadius2D * 1e+4);
            float xiccdcaXY = 1e+10, xiccdcaZ = 1e+10;
            if (xiccTrack.propagateToDCA(primaryVertex, magneticField, &dcaInfo)) {
              xiccdcaXY = dcaInfo.getY();
              xiccdcaZ = dcaInfo.getZ();
            }

            if (std::fabs(xiccdcaXY) > xiccMaxDCAxy || std::fabs(xiccdcaZ) > xiccMaxDCAz) {
              continue; // not pointing to PV
            }

            histos.fill(HIST("hDCAxyXiCC"), std::fabs(xiccdcaXY * 1e+4));
            histos.fill(HIST("hDCAzXiCC"), std::fabs(xiccdcaZ * 1e+4));

            if (std::fabs(thisXiCCcandidate.eta) > xiccMaxEta) {
              continue; // not in central barrel
            }

            histos.fill(HIST("hCharmBuilding"), 3.0f);
            histos.fill(HIST("hMassXiCC"), thisXiCCcandidate.mass);
            GET_HIST(TH1, histPath + "hMassXiCC")->Fill(thisXiCCcandidate.mass);

            histos.fill(HIST("hPtXiCC"), thisXiCCcandidate.pt);
            histos.fill(HIST("hEtaXiCC"), thisXiCCcandidate.eta);
            histos.fill(HIST("h3dMassXiCC"), thisXiCCcandidate.pt, thisXiCCcandidate.eta, thisXiCCcandidate.mass);

            // produce multi-charm table for posterior analysis
            if (fillDerivedTable) {
              multiCharmIdx(
                xiCand.globalIndex(),
                pi1c.globalIndex(), pi2c.globalIndex(),
                picc.globalIndex());

              multiCharmCore(
                thisXiCCcandidate.mass, thisXiCCcandidate.pt,
                thisXiCCcandidate.eta, thisXiCCcandidate.dca,
                thisXiCcandidate.mass, thisXiCcandidate.pt,
                thisXiCcandidate.eta, thisXiCcandidate.dca,
                xi.dcaXY(), xi.dcaZ(),
                xicdcaXY, xicdcaZ,
                xiccdcaXY, xiccdcaZ,
                pi1c.dcaXY(), pi1c.dcaZ(),
                pi2c.dcaXY(), pi2c.dcaZ(),
                picc.dcaXY(), picc.dcaZ(),
                xicDecayRadius2D, xiccDecayRadius2D,
                xicProperLength,
                xicDecayDistanceFromPV,
                xiccProperLength,
                pi1c.pt(), pi2c.pt(), picc.pt(),
                collision.lutConfigId());

              multiCharmPID(
                pi1cTOFDiffInner, pi1c.nSigmaPionInnerTOF(),
                pi1cTOFDiffOuter, pi1c.nSigmaPionOuterTOF(),
                pi1c.hasSigPi(), pi1c.nSigmaPionRich(),
                getPdgCodeForTrack(pi1c),
                pi2cTOFDiffInner, pi2c.nSigmaPionInnerTOF(),
                pi2cTOFDiffOuter, pi2c.nSigmaPionOuterTOF(),
                pi2c.hasSigPi(), pi2c.nSigmaPionRich(),
                getPdgCodeForTrack(pi2c),
                piccTOFDiffInner, picc.nSigmaPionInnerTOF(),
                piccTOFDiffOuter, picc.nSigmaPionOuterTOF(),
                picc.hasSigPi(), picc.nSigmaPionRich(),
                getPdgCodeForTrack(picc));

              multiCharmExtra(
                bach.pt(), bach.eta(),
                bach.dcaXY(), bach.dcaZ(),
                pos.pt(), pos.eta(),
                pos.dcaXY(), pos.dcaZ(),
                neg.pt(), neg.eta(),
                neg.dcaXY(), neg.dcaZ(),
                pi1c.eta(), pi2c.eta(), picc.eta());

              histos.fill(HIST("hPi1cDCAxy"), std::abs(pi1c.dcaXY() * 1e+4));
              histos.fill(HIST("hPi1cDCAz"), std::abs(pi1c.dcaZ() * 1e+4));
              histos.fill(HIST("hPi2cDCAxy"), std::abs(pi2c.dcaXY() * 1e+4));
              histos.fill(HIST("hPi2cDCAz"), std::abs(pi2c.dcaZ() * 1e+4));
              histos.fill(HIST("hPiccDCAxy"), std::abs(picc.dcaXY() * 1e+4));
              histos.fill(HIST("hPiccDCAz"), std::abs(picc.dcaZ() * 1e+4));
            }
          }
          histos.fill(HIST("hCombinationsXiCC"), nCombinationsCC);
        }
      }
      histos.fill(HIST("hCombinationsXiC"), nCombinationsC);
    }
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
  PROCESS_SWITCH(Alice3MulticharmFinder, processGenerated, "fill MC-only histograms", true);
  PROCESS_SWITCH(Alice3MulticharmFinder, processFindXiCC, "find XiCC baryons", true);
  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Alice3MulticharmFinder>(cfgc)};
}
