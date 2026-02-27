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
#include "Common/Core/TableHelper.h"
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
#define INSERT_HIST(name, ...) histPointers[name] = histos.add((name).c_str(), __VA_ARGS__);

using Alice3Tracks = soa::Join<aod::Tracks, aod::TracksCov, aod::Alice3DecayMaps, aod::McTrackLabels, aod::TracksDCA, aod::TracksExtraA3>;
using Alice3Collision = soa::Join<aod::Collisions, aod::OTFLUTConfigId>::iterator;

struct Alice3MulticharmFinder {
  SliceCache cache;

  Produces<aod::MCharmIndices> multiCharmIdx;
  Produces<aod::MCharmCores> multiCharmCore;
  Produces<aod::MCharmExtra> multiCharmExtra;

  // Operation and minimisation criteria

  struct : ConfigurableGroup {
    std::string prefix = "derivedTable"; // JSON group name
    Configurable<bool> fillMCharmIdx{"fillMCharmIdx", true, "fill MCharmIdx[] tables (careful: memory)"};
    Configurable<bool> fillMCharmCore{"fillMCharmCore", true, "fill MCharmCores[] tables (careful: memory)"};
    Configurable<bool> fillMCharmExtra{"fillMCharmExtra", false, "fill MCharmExtra[] tables (careful: memory)"};
  } derivedTable; // allows for gap between peak and bg in case someone wants to

  Configurable<float> cfgMagneticField{"cfgMagneticField", 20.0f, "Magnetic field (in kilogauss) if value not found from geo provider"};
  Configurable<bool> doDCAplots{"doDCAplots", true, "do daughter prong DCA plots for D mesons"};
  Configurable<bool> mcSameMotherCheck{"mcSameMotherCheck", true, "check if tracks come from the same MC mother"};
  Configurable<std::vector<float>> minNTracks{"minNTracks", {-1}, "Minimum number of tracks"};

  Configurable<float> xiMinConstDCAxy{"xiMinConstDCAxy", 0.0005f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinConstDCAz{"xiMinConstDCAz", 0.0005f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinPtDepDCAxy{"xiMinPtDepDCAxy", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinPtDepDCAz{"xiMinPtDepDCAz", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinDecayRadius{"xiMinDecayRadius", 0.5, "Minimum R2D for XiC decay (cm)"};
  Configurable<float> xiMassWindow{"xiMassWindow", 0.005, "Mass window around Xi peak (GeV/c)"};

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
  static constexpr float ToMicrons = 1e+4;

  // filter expressions for pions
  static constexpr uint32_t TrackSelectionPic = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiFromXiC;
  static constexpr uint32_t TrackSelectionPicc = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiFromXiCC;
  float magneticField{};

  // partitions
  Partition<aod::McParticles> trueXi = aod::mcparticle::pdgCode == static_cast<int>(PDG_t::kXiMinus);
  Partition<aod::McParticles> trueXiC = aod::mcparticle::pdgCode == static_cast<int>(o2::constants::physics::kXiCPlus);
  Partition<aod::McParticles> trueXiCC = aod::mcparticle::pdgCode == static_cast<int>(o2::constants::physics::kXiCCPlusPlus);

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
  } thisXicCandidate;

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
  } thisXiccCandidate;

  struct ProngInfo {
    float pt = 1e+10;
    float eta = 1e+10;
    float dcaXY = 1e+10;
    float dcaZ = 1e+10;
  };

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

    const u_int8_t fitterStatusCode = fitter.getFitStatus();
    histos.fill(HIST("hFitterStatusCode"), fitterStatusCode);
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}

    o2::track::TrackParCov t0new = fitter.getTrack(0);
    o2::track::TrackParCov t1new = fitter.getTrack(1);
    t0new.getPxPyPzGlo(thisXiccCandidate.prong0mom);
    t1new.getPxPyPzGlo(thisXiccCandidate.prong1mom);

    // get decay vertex coordinates
    const auto& vtx = fitter.getPCACandidate();
    for (size_t i = 0; i < thisXiccCandidate.xyz.size(); i++) {
      thisXiccCandidate.xyz[i] = vtx[i];
    }

    // compute cov mat
    for (int ii = 0; ii < o2::track::kLabCovMatSize; ii++) {
      thisXiccCandidate.parentTrackCovMatrix[ii] = 0.0f;
    }

    std::array<float, o2::track::kLabCovMatSize> covA = {0};
    std::array<float, o2::track::kLabCovMatSize> covB = {0};
    fitter.getTrack(0).getCovXYZPxPyPzGlo(covA);
    fitter.getTrack(1).getCovXYZPxPyPzGlo(covB);

    for (size_t i = 0; i < MomentumIndices.size(); i++) {
      int j = MomentumIndices[i];
      thisXiccCandidate.parentTrackCovMatrix[j] = covA[j] + covB[j];
    }

    auto covVtx = fitter.calcPCACovMatrix();
    thisXiccCandidate.parentTrackCovMatrix[0] = covVtx(0, 0);
    thisXiccCandidate.parentTrackCovMatrix[1] = covVtx(1, 0);
    thisXiccCandidate.parentTrackCovMatrix[2] = covVtx(1, 1);
    thisXiccCandidate.parentTrackCovMatrix[3] = covVtx(2, 0);
    thisXiccCandidate.parentTrackCovMatrix[4] = covVtx(2, 1);
    thisXiccCandidate.parentTrackCovMatrix[5] = covVtx(2, 2);

    // set relevant values
    thisXiccCandidate.dca = std::sqrt(fitter.getChi2AtPCACandidate());
    if (thisXiccCandidate.dca > xiccMaxDauDCA) {
      return false;
    }

    thisXiccCandidate.mass = RecoDecay::m(std::array{std::array{thisXiccCandidate.prong0mom[0], thisXiccCandidate.prong0mom[1], thisXiccCandidate.prong0mom[2]}, std::array{thisXiccCandidate.prong1mom[0], thisXiccCandidate.prong1mom[1], thisXiccCandidate.prong1mom[2]}}, std::array{mass0, mass1});

    if (std::fabs(thisXiccCandidate.mass - o2::constants::physics::MassXiCCPlusPlus) > xiccMassWindow) {
      return false;
    }

    thisXiccCandidate.pt = std::hypot(thisXiccCandidate.prong0mom[0] + thisXiccCandidate.prong1mom[0], thisXiccCandidate.prong0mom[1] + thisXiccCandidate.prong1mom[1]);
    thisXiccCandidate.eta = RecoDecay::eta(std::array{thisXiccCandidate.prong0mom[0] + thisXiccCandidate.prong1mom[0], thisXiccCandidate.prong0mom[1] + thisXiccCandidate.prong1mom[1], thisXiccCandidate.prong0mom[2] + thisXiccCandidate.prong1mom[2]});
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

    const u_int8_t fitter3StatusCode = fitter3.getFitStatus();
    histos.fill(HIST("hFitter3StatusCode"), fitter3StatusCode);
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}

    t0 = fitter3.getTrack(0);
    t1 = fitter3.getTrack(1);
    t2 = fitter3.getTrack(2);
    t0.getPxPyPzGlo(thisXicCandidate.prong0mom);
    t1.getPxPyPzGlo(thisXicCandidate.prong1mom);
    t2.getPxPyPzGlo(thisXicCandidate.prong2mom);

    // get decay vertex coordinates
    const auto& vtx = fitter3.getPCACandidate();
    for (size_t i = 0; i < thisXicCandidate.xyz.size(); i++) {
      thisXicCandidate.xyz[i] = vtx[i];
    }

    // compute cov mat
    for (int ii = 0; ii < o2::track::kLabCovMatSize; ii++) {
      thisXicCandidate.parentTrackCovMatrix[ii] = 0.0f;
    }

    std::array<float, o2::track::kLabCovMatSize> covA = {0};
    std::array<float, o2::track::kLabCovMatSize> covB = {0};
    std::array<float, o2::track::kLabCovMatSize> covC = {0};
    fitter3.getTrack(0).getCovXYZPxPyPzGlo(covA);
    fitter3.getTrack(1).getCovXYZPxPyPzGlo(covB);
    fitter3.getTrack(2).getCovXYZPxPyPzGlo(covC);

    for (size_t i = 0; i < MomentumIndices.size(); i++) {
      int j = MomentumIndices[i];
      thisXicCandidate.parentTrackCovMatrix[j] = covA[j] + covB[j] + covC[j];
    }

    auto covVtx = fitter3.calcPCACovMatrix();
    thisXicCandidate.parentTrackCovMatrix[0] = covVtx(0, 0);
    thisXicCandidate.parentTrackCovMatrix[1] = covVtx(1, 0);
    thisXicCandidate.parentTrackCovMatrix[2] = covVtx(1, 1);
    thisXicCandidate.parentTrackCovMatrix[3] = covVtx(2, 0);
    thisXicCandidate.parentTrackCovMatrix[4] = covVtx(2, 1);
    thisXicCandidate.parentTrackCovMatrix[5] = covVtx(2, 2);

    // set relevant values
    thisXicCandidate.dca = std::sqrt(fitter3.getChi2AtPCACandidate());
    if (thisXicCandidate.dca > xicMaxDauDCA) {
      return false;
    }
    thisXicCandidate.mass = RecoDecay::m(std::array{std::array{thisXicCandidate.prong0mom[0], thisXicCandidate.prong0mom[1], thisXicCandidate.prong0mom[2]}, std::array{thisXicCandidate.prong1mom[0], thisXicCandidate.prong1mom[1], thisXicCandidate.prong1mom[2]}, std::array{thisXicCandidate.prong2mom[0], thisXicCandidate.prong2mom[1], thisXicCandidate.prong2mom[2]}}, std::array{p0mass, p1mass, p2mass});
    thisXicCandidate.pt = std::hypot(thisXicCandidate.prong0mom[0] + thisXicCandidate.prong1mom[0] + thisXicCandidate.prong2mom[0], thisXicCandidate.prong0mom[1] + thisXicCandidate.prong1mom[1] + thisXicCandidate.prong2mom[1]);
    thisXicCandidate.eta = RecoDecay::eta(std::array{thisXicCandidate.prong0mom[0] + thisXicCandidate.prong1mom[0] + thisXicCandidate.prong2mom[0], thisXicCandidate.prong0mom[1] + thisXicCandidate.prong1mom[1] + thisXicCandidate.prong2mom[1], thisXicCandidate.prong0mom[2] + thisXicCandidate.prong1mom[2] + thisXicCandidate.prong2mom[2]});
    return true;
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

  void init(o2::framework::InitContext& initContext)
  {
    const bool foundMagneticField = common::core::getTaskOptionValue(initContext, "on-the-fly-detector-geometry-provider", "magneticField", magneticField, false);
    if (!foundMagneticField) {
      LOG(info) << "Could not retrieve magnetic field from geometry provider.";
      LOG(info) << "Using value from configurable cfgMagneticField: " << cfgMagneticField;
      magneticField = cfgMagneticField;
    } else {
      LOG(info) << "Using magnetic field form geometry provider with value: " << magneticField;
    }

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

    auto hFitterStatusCode = histos.add<TH1>("hFitterStatusCode", "hFitterStatusCode", kTH1D, {{15, -0.5, 14.5}});
    hFitterStatusCode->GetXaxis()->SetBinLabel(1, "None"); // no status set (should not be possible!)

    /* Good Conditions */
    hFitterStatusCode->GetXaxis()->SetBinLabel(2, "Converged"); // fit converged
    hFitterStatusCode->GetXaxis()->SetBinLabel(3, "MaxIter");   // max iterations reached before fit convergence

    /* Error Conditions */
    hFitterStatusCode->GetXaxis()->SetBinLabel(4, "NoCrossing");       // no reasaonable crossing was found
    hFitterStatusCode->GetXaxis()->SetBinLabel(5, "RejRadius");        // radius of crossing was not acceptable
    hFitterStatusCode->GetXaxis()->SetBinLabel(6, "RejTrackX");        // one candidate track x was below the mimimum required radius
    hFitterStatusCode->GetXaxis()->SetBinLabel(7, "RejTrackRoughZ");   // rejected by rough cut on tracks Z difference
    hFitterStatusCode->GetXaxis()->SetBinLabel(8, "RejChi2Max");       // rejected by maximum chi2 cut
    hFitterStatusCode->GetXaxis()->SetBinLabel(9, "FailProp");         // propagation of at least prong to PCA failed
    hFitterStatusCode->GetXaxis()->SetBinLabel(10, "FailInvCov");      // inversion of cov.-matrix failed
    hFitterStatusCode->GetXaxis()->SetBinLabel(11, "FailInvWeight");   // inversion of Ti weight matrix failed
    hFitterStatusCode->GetXaxis()->SetBinLabel(12, "FailInv2ndDeriv"); // inversion of 2nd derivatives failed
    hFitterStatusCode->GetXaxis()->SetBinLabel(13, "FailCorrTracks");  // correction of tracks to updated x failed
    hFitterStatusCode->GetXaxis()->SetBinLabel(14, "FailCloserAlt");   // alternative PCA is closer
    hFitterStatusCode->GetXaxis()->SetBinLabel(15, "NStatusesDefined");

    auto hFitter3StatusCode = histos.add<TH1>("hFitter3StatusCode", "hFitter3StatusCode", kTH1D, {{15, -0.5, 14.5}});
    hFitter3StatusCode->GetXaxis()->SetBinLabel(1, "None"); // no status set (should not be possible!)

    /* Good Conditions */
    hFitter3StatusCode->GetXaxis()->SetBinLabel(2, "Converged"); // fit converged
    hFitter3StatusCode->GetXaxis()->SetBinLabel(3, "MaxIter");   // max iterations reached before fit convergence

    /* Error Conditions */
    hFitter3StatusCode->GetXaxis()->SetBinLabel(4, "NoCrossing");       // no reasaonable crossing was found
    hFitter3StatusCode->GetXaxis()->SetBinLabel(5, "RejRadius");        // radius of crossing was not acceptable
    hFitter3StatusCode->GetXaxis()->SetBinLabel(6, "RejTrackX");        // one candidate track x was below the mimimum required radius
    hFitter3StatusCode->GetXaxis()->SetBinLabel(7, "RejTrackRoughZ");   // rejected by rough cut on tracks Z difference
    hFitter3StatusCode->GetXaxis()->SetBinLabel(8, "RejChi2Max");       // rejected by maximum chi2 cut
    hFitter3StatusCode->GetXaxis()->SetBinLabel(9, "FailProp");         // propagation of at least prong to PCA failed
    hFitter3StatusCode->GetXaxis()->SetBinLabel(10, "FailInvCov");      // inversion of cov.-matrix failed
    hFitter3StatusCode->GetXaxis()->SetBinLabel(11, "FailInvWeight");   // inversion of Ti weight matrix failed
    hFitter3StatusCode->GetXaxis()->SetBinLabel(12, "FailInv2ndDeriv"); // inversion of 2nd derivatives failed
    hFitter3StatusCode->GetXaxis()->SetBinLabel(13, "FailCorrTracks");  // correction of tracks to updated x failed
    hFitter3StatusCode->GetXaxis()->SetBinLabel(14, "FailCloserAlt");   // alternative PCA is closer
    hFitter3StatusCode->GetXaxis()->SetBinLabel(15, "NStatusesDefined");

    INSERT_HIST(std::string("h2dGenXi"), "h2dGenXi", {kTH2D, {{axisPt, axisEta}}});
    INSERT_HIST(std::string("h2dGenXiC"), "h2dGenXiC", {kTH2D, {{axisPt, axisEta}}});
    INSERT_HIST(std::string("h2dGenXiCC"), "h2dGenXiCC", {kTH2D, {{axisPt, axisEta}}});
  }

  void initDetectorConfiguration(const int icfg)
  {
    if (std::find(savedConfigs.begin(), savedConfigs.end(), icfg) != savedConfigs.end()) {
      return;
    }

    savedConfigs.push_back(icfg);
    const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";

    // This histogram bookkeeps the attempts at DCA minimization and their eventual
    // failure rates.
    // --- 0: attempt XiC, 1: success XiC
    // --- 2: attempt XiCC, 3: success XiCC
    INSERT_HIST(histPath + "hCharmBuilding", "hCharmBuilding", {kTH1D, {{10, -0.5, 9.5f}}});
    INSERT_HIST(histPath + "hMultiCharmBuilding", "hMultiCharmBuilding", {kTH1D, {{10, -0.5, 9.5f}}});

    INSERT_HIST(histPath + "hMassXi", "hMassXi", {kTH1D, {{axisXiMass}}});
    INSERT_HIST(histPath + "hMassXiC", "hMassXiC", {kTH1D, {{axisXiCMass}}});

    INSERT_HIST(histPath + "hEtaXiCC", "hEtaXiCC", {kTH1D, {{axisEta}}});
    INSERT_HIST(histPath + "hPtXiCC", "hPtXiCC", {kTH1D, {{axisPt}}});
    INSERT_HIST(histPath + "h3dMassXiCC", "h3dMassXiCC", {kTH3D, {{axisPt, axisEta, axisXiCCMass}}});

    INSERT_HIST(histPath + "hDCAXiCDaughters", "hDCAXiCDaughters", {kTH1D, {{axisDCAXiCDaughters}}});
    INSERT_HIST(histPath + "hDCAXiCCDaughters", "hDCAXiCCDaughters", {kTH1D, {{axisDCAXiCCDaughters}}});
    INSERT_HIST(histPath + "hDCAxyXi", "hDCAxyXi", {kTH1D, {{axisDCA}}});
    INSERT_HIST(histPath + "hDCAzXi", "hDCAzXi", {kTH1D, {{axisDCA}}});

    INSERT_HIST(histPath + "hDCAxyXiC", "hDCAxyXiC", {kTH1D, {{axisDCA}}});
    INSERT_HIST(histPath + "hDCAzXiC", "hDCAzXiC", {kTH1D, {{axisDCA}}});

    INSERT_HIST(histPath + "hDCAxyXiCC", "hDCAxyXiCC", {kTH1D, {{axisDCA}}});
    INSERT_HIST(histPath + "hDCAzXiCC", "hDCAzXiCC", {kTH1D, {{axisDCA}}});

    INSERT_HIST(histPath + "hPi1cPt", "hPi1cPt", {kTH1D, {{axisPt}}});
    INSERT_HIST(histPath + "hPi2cPt", "hPi2cPt", {kTH1D, {{axisPt}}});
    INSERT_HIST(histPath + "hPiccPt", "hPiccPt", {kTH1D, {{axisPt}}});

    INSERT_HIST(histPath + "hPi1cDCAxy", "hPi1cDCAxy", {kTH1D, {{axisDCA}}});
    INSERT_HIST(histPath + "hPi1cDCAz", "hPi1cDCAz", {kTH1D, {{axisDCA}}});
    INSERT_HIST(histPath + "hPi2cDCAxy", "hPi2cDCAxy", {kTH1D, {{axisDCA}}});
    INSERT_HIST(histPath + "hPi2cDCAz", "hPi2cDCAz", {kTH1D, {{axisDCA}}});
    INSERT_HIST(histPath + "hPiccDCAxy", "hPiccDCAxy", {kTH1D, {{axisDCA}}});
    INSERT_HIST(histPath + "hPiccDCAz", "hPiccDCAz", {kTH1D, {{axisDCA}}});

    INSERT_HIST(histPath + "hMinXiDecayRadius", "hMinXiDecayRadius", {kTH1D, {{axisRadius2DXi}}});
    INSERT_HIST(histPath + "hMinXiCDecayRadius", "hMinXiCDecayRadius", {kTH1D, {{axisRadius}}});
    INSERT_HIST(histPath + "hMinXiCCDecayRadius", "hMinXiCCDecayRadius", {kTH1D, {{axisRadius}}});

    INSERT_HIST(histPath + "hMinXicDecayDistanceFromPV", "hMinXicDecayDistanceFromPV", {kTH1D, {{axisDecayLength}}});
    INSERT_HIST(histPath + "hProperLengthXiC", "hProperLengthXiC", {kTH1D, {{axisDecayLength}}});
    INSERT_HIST(histPath + "hProperLengthXiCC", "hProperLengthXiCC", {kTH1D, {{axisDecayLength}}});

    INSERT_HIST(histPath + "hInnerTOFTrackTimeRecoPi1c", "hInnerTOFTrackTimeRecoPi1c", {kTH1D, {{axisTOFTrack}}});
    INSERT_HIST(histPath + "hInnerTOFTrackTimeRecoPi2c", "hInnerTOFTrackTimeRecoPi2c", {kTH1D, {{axisTOFTrack}}});
    INSERT_HIST(histPath + "hInnerTOFTrackTimeRecoPicc", "hInnerTOFTrackTimeRecoPicc", {kTH1D, {{axisTOFTrack}}});

    INSERT_HIST(histPath + "hOuterTOFTrackTimeRecoPi1c", "hOuterTOFTrackTimeRecoPi1c", {kTH1D, {{axisTOFTrack}}});
    INSERT_HIST(histPath + "hOuterTOFTrackTimeRecoPi2c", "hOuterTOFTrackTimeRecoPi2c", {kTH1D, {{axisTOFTrack}}});
    INSERT_HIST(histPath + "hOuterTOFTrackTimeRecoPicc", "hOuterTOFTrackTimeRecoPicc", {kTH1D, {{axisTOFTrack}}});

    INSERT_HIST(histPath + "hXiRadiusVsXicRadius", "hXiRadiusVsXicRadius", {kTH2D, {{axisRadius2D, axisRadius2D}}});
    INSERT_HIST(histPath + "hXicRadiusVsXiccRadius", "hXicRadiusVsXiccRadius", {kTH2D, {{axisRadius2D, axisRadius2D}}});

    INSERT_HIST(histPath + "hMassXiCC", "hMassXiCC", {kTH1D, {{axisXiCCMass}}});
    INSERT_HIST(histPath + "hNCollisions", "hNCollisions", {kTH1D, {{2, 0.5, 2.5}}});
    INSERT_HIST(histPath + "hNTracks", "hNTracks", {kTH1D, {{20000, 0, 20000}}});

    // These histograms bookkeep the exact number of combinations attempted
    // CombinationsXiC: triplets Xi-pi-pi considered per Xi
    // CombinationsXiCC: doublets XiC-pi considered per XiC
    INSERT_HIST(histPath + "hCombinationsXiC", "hCombinationsXiC", {kTH1D, {{axisNConsidered}}});
    INSERT_HIST(histPath + "hCombinationsXiCC", "hCombinationsXiCC", {kTH1D, {{axisNConsidered}}});

    if (doDCAplots) {
      INSERT_HIST(histPath + "h2dDCAxyVsPtXiFromXiC", "h2dDCAxyVsPtXiFromXiC", {kTH2D, {{axisPt, axisDCA2D}}});
      INSERT_HIST(histPath + "h2dDCAxyVsPtPiFromXiC", "h2dDCAxyVsPtPiFromXiC", {kTH2D, {{axisPt, axisDCA2D}}});
      INSERT_HIST(histPath + "h2dDCAxyVsPtPiFromXiCC", "h2dDCAxyVsPtPiFromXiCC", {kTH2D, {{axisPt, axisDCA2D}}});

      INSERT_HIST(histPath + "h2dDCAzVsPtXiFromXiC", "h2dDCAzVsPtXiFromXiC", {kTH2D, {{axisPt, axisDCA2D}}});
      INSERT_HIST(histPath + "h2dDCAzVsPtPiFromXiC", "h2dDCAzVsPtPiFromXiC", {kTH2D, {{axisPt, axisDCA2D}}});
      INSERT_HIST(histPath + "h2dDCAzVsPtPiFromXiCC", "h2dDCAzVsPtPiFromXiCC", {kTH2D, {{axisPt, axisDCA2D}}});
    }
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processGenerated(aod::McParticles const&)
  {
    for (auto const& mcParticle : trueXi) {
      GET_HIST(TH2, "h2dGenXi")->Fill(mcParticle.pt(), mcParticle.eta());
    }
    for (auto const& mcParticle : trueXiC) {
      GET_HIST(TH2, "h2dGenXiC")->Fill(mcParticle.pt(), mcParticle.eta());
    }
    for (auto const& mcParticle : trueXiCC) {
      GET_HIST(TH2, "h2dGenXiCC")->Fill(mcParticle.pt(), mcParticle.eta());
    }
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFindXicc(Alice3Collision const& collision, Alice3Tracks const& tracks, aod::McParticles const&, aod::UpgradeCascades const& cascades)
  {
    const std::string histPath = "Configuration_" + std::to_string(collision.lutConfigId()) + "/";
    initDetectorConfiguration(collision.lutConfigId());

    GET_HIST(TH1, histPath + "hNCollisions")->Fill(1);
    GET_HIST(TH1, histPath + "hNTracks")->Fill(tracks.size());
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
        GET_HIST(TH2, histPath + "h2dDCAxyVsPtPiFromXiC")->Fill(track.pt(), track.dcaXY() * ToMicrons);
        GET_HIST(TH2, histPath + "h2dDCAxyVsPtPiFromXiC")->Fill(track.pt(), track.dcaZ() * ToMicrons);
      }
      if (BIT_CHECK(track.decayMap(), kTruePiFromXiCC)) {
        GET_HIST(TH2, histPath + "h2dDCAxyVsPtPiFromXiCC")->Fill(track.pt(), track.dcaXY() * ToMicrons);
        GET_HIST(TH2, histPath + "h2dDCAxyVsPtPiFromXiCC")->Fill(track.pt(), track.dcaZ() * ToMicrons);
      }
    }

    for (auto const& xiCand : cascades) {
      auto xi = xiCand.cascadeTrack_as<Alice3Tracks>(); // de-reference cascade track
      GET_HIST(TH1, histPath + "hMassXi")->Fill(xiCand.mXi());
      GET_HIST(TH2, histPath + "h2dDCAxyVsPtXiFromXiC")->Fill(xi.pt(), xi.dcaXY() * ToMicrons);
      GET_HIST(TH2, histPath + "h2dDCAzVsPtXiFromXiC")->Fill(xi.pt(), xi.dcaZ() * ToMicrons);

      if (std::fabs(xiCand.mXi() - o2::constants::physics::MassXiMinus) > xiMassWindow) {
        continue; // out of mass region
      }

      uint32_t nCombinationsC = 0;
      if (!BIT_CHECK(xi.decayMap(), kTrueXiFromXiC)) {
        continue;
      }

      if (std::fabs(xi.dcaXY()) < xiMinConstDCAxy || std::fabs(xi.dcaZ()) < xiMinConstDCAz) {
        continue; // likely a primary xi
      }

      GET_HIST(TH1, histPath + "hDCAxyXi")->Fill(xi.dcaXY() * ToMicrons);
      GET_HIST(TH1, histPath + "hDCAzXi")->Fill(xi.dcaZ() * ToMicrons);
      if (xiCand.cascRadius() < xiMinDecayRadius) {
        continue;
      }

      GET_HIST(TH1, histPath + "hMinXiDecayRadius")->Fill(xiCand.cascRadius());
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

        GET_HIST(TH1, histPath + "hPi1cPt")->Fill(pi1c.pt());
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

          GET_HIST(TH1, histPath + "hPi2cPt")->Fill(pi2c.pt());
          // if I am here, it means this is a triplet to be considered for XiC vertexing.
          // will now attempt to build a three-body decay candidate with these three track rows.
          nCombinationsC++;
          GET_HIST(TH1, histPath + "hCharmBuilding")->Fill(0.0f);

          if (!buildDecayCandidateThreeBody(xi, pi1c, pi2c, o2::constants::physics::MassXiMinus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged)) {
            continue; // failed at building candidate
          }

          GET_HIST(TH1, histPath + "hDCAXiCDaughters")->Fill(thisXicCandidate.dca * ToMicrons);
          if (std::fabs(thisXicCandidate.mass - o2::constants::physics::MassXiCPlus) > xicMassWindow) {
            continue; // out of mass region
          }

          GET_HIST(TH1, histPath + "hCharmBuilding")->Fill(1.0f);
          const std::array<float, 3> momentumC = {
            thisXicCandidate.prong0mom[0] + thisXicCandidate.prong1mom[0] + thisXicCandidate.prong2mom[0],
            thisXicCandidate.prong0mom[1] + thisXicCandidate.prong1mom[1] + thisXicCandidate.prong2mom[1],
            thisXicCandidate.prong0mom[2] + thisXicCandidate.prong1mom[2] + thisXicCandidate.prong2mom[2]};

          o2::track::TrackParCov xicTrack(thisXicCandidate.xyz, momentumC, thisXicCandidate.parentTrackCovMatrix, +1);
          float xicDecayRadius2D = std::hypot(thisXicCandidate.xyz[0], thisXicCandidate.xyz[1]);
          if (xicDecayRadius2D < xiccMinDecayRadius) {
            continue; // do not take if radius too small, likely a primary combination
          }

          GET_HIST(TH1, histPath + "hCharmBuilding")->Fill(2.0f);
          GET_HIST(TH1, histPath + "hMinXiCDecayRadius")->Fill(xicDecayRadius2D * ToMicrons);
          if (xicDecayRadius2D > xiCand.cascRadius()) {
            continue;
          }

          GET_HIST(TH1, histPath + "hCharmBuilding")->Fill(3.0f);
          GET_HIST(TH2, histPath + "hXiRadiusVsXicRadius")->Fill(xiCand.cascRadius() * ToMicrons, xicDecayRadius2D * ToMicrons);
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

          GET_HIST(TH1, histPath + "hCharmBuilding")->Fill(4.0f);
          GET_HIST(TH1, histPath + "hDCAxyXiC")->Fill(std::fabs(xicdcaXY * ToMicrons));
          GET_HIST(TH1, histPath + "hDCAzXiC")->Fill(std::fabs(xicdcaZ * ToMicrons));
          GET_HIST(TH1, histPath + "hMassXiC")->Fill(thisXicCandidate.mass);

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

            GET_HIST(TH1, histPath + "hMultiCharmBuilding")->Fill(0.0f);
            GET_HIST(TH1, histPath + "hPiccPt")->Fill(picc.pt());
            o2::track::TrackParCov piccTrack = getTrackParCov(picc);
            nCombinationsCC++;
            if (!buildDecayCandidateTwoBody(xicTrack, piccTrack, o2::constants::physics::MassXiCPlus, o2::constants::physics::MassPionCharged)) {
              continue; // failed at building candidate
            }

            GET_HIST(TH1, histPath + "hMultiCharmBuilding")->Fill(1.0f);
            GET_HIST(TH1, histPath + "hDCAXiCCDaughters")->Fill(thisXiccCandidate.dca * ToMicrons);
            const std::array<float, 3> momentumCC = {
              thisXiccCandidate.prong0mom[0] + thisXiccCandidate.prong1mom[0],
              thisXiccCandidate.prong0mom[1] + thisXiccCandidate.prong1mom[1],
              thisXiccCandidate.prong0mom[2] + thisXiccCandidate.prong1mom[2]};

            o2::track::TrackParCov xiccTrack(thisXiccCandidate.xyz, momentumCC, thisXiccCandidate.parentTrackCovMatrix, +2);
            float xiccDecayRadius2D = std::hypot(thisXiccCandidate.xyz[0], thisXiccCandidate.xyz[1]);
            if (xiccDecayRadius2D < xiccMinDecayRadius) {
              continue; // do not take if radius too small, likely a primary combination
            }

            GET_HIST(TH1, histPath + "hMultiCharmBuilding")->Fill(2.0f);
            GET_HIST(TH1, histPath + "hDCAXiCCDaughters")->Fill(thisXiccCandidate.dca * ToMicrons);
            float totalMomentumC = std::hypot(momentumC[0], momentumC[1], momentumC[2]);
            float decayLengthXiC = std::hypot(
              thisXicCandidate.xyz[0] - thisXiccCandidate.xyz[0],
              thisXicCandidate.xyz[1] - thisXiccCandidate.xyz[1],
              thisXicCandidate.xyz[2] - thisXiccCandidate.xyz[2]);
            float xicProperLength = decayLengthXiC * thisXicCandidate.mass / totalMomentumC;

            if (xicProperLength < xicMinProperLength || xicProperLength > xicMaxProperLength) {
              continue; // likely background
            }

            GET_HIST(TH1, histPath + "hMultiCharmBuilding")->Fill(3.0f);
            GET_HIST(TH1, histPath + "hProperLengthXiC")->Fill(xicProperLength * ToMicrons);
            float xicDistanceFromPV = std::hypot(
              thisXicCandidate.xyz[0] - collision.posX(),
              thisXicCandidate.xyz[1] - collision.posY(),
              thisXicCandidate.xyz[2] - collision.posZ());
            float xicDecayDistanceFromPV = xicDistanceFromPV * thisXicCandidate.mass / totalMomentumC;
            if (xicDecayDistanceFromPV < xicMinDecayDistanceFromPV) {
              continue; // too close to PV
            }

            GET_HIST(TH1, histPath + "hMultiCharmBuilding")->Fill(4.0f);
            GET_HIST(TH1, histPath + "hMinXicDecayDistanceFromPV")->Fill(xicDecayDistanceFromPV * ToMicrons);
            float totalMomentumCC = std::hypot(momentumCC[0], momentumCC[1], momentumCC[2]);
            float decayLengthXiCC = std::hypot(
              thisXiccCandidate.xyz[0] - collision.posX(),
              thisXiccCandidate.xyz[1] - collision.posY(),
              thisXiccCandidate.xyz[2] - collision.posZ());
            float xiccProperLength = decayLengthXiCC * thisXiccCandidate.mass / totalMomentumCC;
            if (xiccProperLength < xiccMinProperLength || xiccProperLength > xicMaxProperLength) {
              continue; // likely background
            }

            GET_HIST(TH1, histPath + "hMultiCharmBuilding")->Fill(5.0f);
            GET_HIST(TH1, histPath + "hProperLengthXiCC")->Fill(xiccProperLength * ToMicrons);
            if (xiccDecayRadius2D > xicDecayRadius2D) {
              continue; // XiCC should decay before XiC
            }

            GET_HIST(TH1, histPath + "hMultiCharmBuilding")->Fill(6.0f);
            GET_HIST(TH2, histPath + "hXicRadiusVsXiccRadius")->Fill(xicDecayRadius2D * ToMicrons, xiccDecayRadius2D * ToMicrons);
            float xiccdcaXY = 1e+10, xiccdcaZ = 1e+10;
            if (xiccTrack.propagateToDCA(primaryVertex, magneticField, &dcaInfo)) {
              xiccdcaXY = dcaInfo.getY();
              xiccdcaZ = dcaInfo.getZ();
            }

            if (std::fabs(xiccdcaXY) > xiccMaxDCAxy || std::fabs(xiccdcaZ) > xiccMaxDCAz) {
              continue; // not pointing to PV
            }

            GET_HIST(TH1, histPath + "hMultiCharmBuilding")->Fill(7.0f);
            GET_HIST(TH1, histPath + "hDCAxyXiCC")->Fill(xiccdcaXY * ToMicrons);
            GET_HIST(TH1, histPath + "hDCAzXiCC")->Fill(xiccdcaZ * ToMicrons);
            if (std::fabs(thisXiccCandidate.eta) > xiccMaxEta) {
              continue; // not in central barrel
            }

            GET_HIST(TH1, histPath + "hMultiCharmBuilding")->Fill(8.0f);
            GET_HIST(TH1, histPath + "hMassXiCC")->Fill(thisXiccCandidate.mass);
            GET_HIST(TH1, histPath + "hPtXiCC")->Fill(thisXiccCandidate.pt);
            GET_HIST(TH1, histPath + "hEtaXiCC")->Fill(thisXiccCandidate.eta);
            GET_HIST(TH3, histPath + "h3dMassXiCC")->Fill(thisXiccCandidate.pt, thisXiccCandidate.eta, thisXiccCandidate.mass);

            GET_HIST(TH1, histPath + "hPi1cDCAxy")->Fill(std::abs(pi1c.dcaXY() * ToMicrons));
            GET_HIST(TH1, histPath + "hPi1cDCAz")->Fill(std::abs(pi1c.dcaZ() * ToMicrons));
            GET_HIST(TH1, histPath + "hPi2cDCAxy")->Fill(std::abs(pi2c.dcaXY() * ToMicrons));
            GET_HIST(TH1, histPath + "hPi2cDCAz")->Fill(std::abs(pi2c.dcaZ() * ToMicrons));
            GET_HIST(TH1, histPath + "hPiccDCAxy")->Fill(std::abs(picc.dcaXY() * ToMicrons));
            GET_HIST(TH1, histPath + "hPiccDCAz")->Fill(std::abs(picc.dcaZ() * ToMicrons));

            // produce multi-charm table for posterior analysis
            if (derivedTable.fillMCharmIdx) {
              multiCharmIdx(
                xiCand.globalIndex(),
                pi1c.globalIndex(), pi2c.globalIndex(),
                picc.globalIndex());
            }

            if (derivedTable.fillMCharmCore) {
              multiCharmCore(
                thisXiccCandidate.mass, thisXiccCandidate.pt,
                thisXiccCandidate.eta, thisXiccCandidate.dca,
                thisXicCandidate.mass, thisXicCandidate.pt,
                thisXicCandidate.eta, thisXicCandidate.dca,
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
            }

            if (derivedTable.fillMCharmExtra) {
              ProngInfo bachelor, positive, negative;
              if (xiCand.has_bachTrack()) {
                auto bach = xiCand.bachTrack_as<Alice3Tracks>(); // de-reference bach track
                bachelor.pt = bach.pt();
                bachelor.eta = bach.eta();
                bachelor.dcaXY = bach.dcaXY();
                bachelor.dcaZ = bach.dcaZ();
              }

              if (xiCand.has_negTrack()) {
                auto neg = xiCand.negTrack_as<Alice3Tracks>(); // de-reference neg track
                negative.pt = neg.pt();
                negative.eta = neg.eta();
                negative.dcaXY = neg.dcaXY();
                negative.dcaZ = neg.dcaZ();
              }

              if (xiCand.has_posTrack()) {
                auto pos = xiCand.posTrack_as<Alice3Tracks>(); // de-reference pos track
                positive.pt = pos.pt();
                positive.eta = pos.eta();
                positive.dcaXY = pos.dcaXY();
                positive.dcaZ = pos.dcaZ();
              }

              multiCharmExtra(
                bachelor.pt, bachelor.eta,
                bachelor.dcaXY, bachelor.dcaZ,
                positive.pt, positive.eta,
                positive.dcaXY, positive.dcaZ,
                negative.pt, negative.eta,
                negative.dcaXY, negative.dcaZ,
                pi1c.eta(), pi2c.eta(), picc.eta());
            }
          }
          GET_HIST(TH1, histPath + "hCombinationsXiCC")->Fill(nCombinationsCC);
        }
      }
      GET_HIST(TH1, histPath + "hCombinationsXiC")->Fill(nCombinationsC);
    }
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
  PROCESS_SWITCH(Alice3MulticharmFinder, processGenerated, "fill MC-only histograms", true);
  PROCESS_SWITCH(Alice3MulticharmFinder, processFindXicc, "find XiCC baryons", true);
  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Alice3MulticharmFinder>(cfgc)};
}
